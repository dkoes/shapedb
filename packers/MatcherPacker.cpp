/*
 * MatcherPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "MatcherPacker.h"
#include <boost/unordered_set.hpp>
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/graph_to_eps.h>
#include <lemon/lp.h>
#include <ANN/ANN.h>
#include "Timer.h"
using namespace lemon;

void MatcherPacker::pack(const DataViewer* dv, vector<Cluster>& clusters) const
{
	if (dv->size() == 0)
		return;

	vector<unsigned> indices(dv->size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
		indices[i] = i;

	//start with each object in its own cluster
	clusters.clear();
	clusters.resize(indices.size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	//combine everything as much as possible
	unsigned curSz = 1;

	if (K == 0)
	{
		FullCache dcache(dv);

		while (curSz < packSize
				&& fullMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;
	}
	else //use knn
	{
		OnDemandCache dcache(dv);

		if (doQuadPack)
		{
			while (curSz < packSize
					&& knnQuadMergeClusters(dv, clusters, curSz, dcache))
				curSz *= 4;

			//always run a cleanup of any small guys
			if (knnMergeClusters(dv, clusters, curSz, dcache))
				curSz *= 2;
		}

		while (curSz < packSize && knnMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;

		//cout << "dcache " << dcache.size() << "\n";
	}

}


/* find the pseudo-optimal matching of every pair of clusters
 * use the knn graph instead of the full graph
 * maxSz is what the largest current cluster should be
 */
bool MatcherPacker::knnMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache) const
{
	unsigned N = clusters.size();
	SmartGraph SG;
	SmartGraph::EdgeMap<double> Sweights(SG);
	SmartGraph::NodeMap<unsigned> Sindex(SG);
	vector<SmartGraph::Node> Snodes;

	ClusterCache ccache(clusters.size());
	makeKNNGraph(D, clusters, maxSz, dcache, ccache, SG, Sweights, Sindex,
			Snodes);

	Timer t;
	MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<double> > matcher(SG,
			Sweights);
	matcher.run();
	//merge clusters that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for (unsigned i = 0; i < N; i++)
	{
		if (!packed[i])
		{
			SmartGraph::Node m = matcher.mate(Snodes[i]);
			if (!SG.valid(m))
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
			else
			{
				unsigned j = Sindex[m];

				if (clusters[i].size() + clusters[j].size() <= packSize) //merge
				{
					newclusters.push_back(Cluster());
					newclusters.back().mergeInto(clusters[i], clusters[j]);
					didmerge = true;
					packed[i] = true;
					packed[j] = true;
				}
				else
				{
					//didn't merge, keep clusters as is
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[i]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[j]);
					packed[i] = true;
					packed[j] = true;
				}
			}
		}
	}

	swap(clusters, newclusters);
	return didmerge;
}

//write dot file from graph
static void writeDot(SmartGraph& g, SmartGraph::EdgeMap<double>& length,
		SmartGraph::NodeMap<unsigned>& index,
		ostream& out)
{
	out << "graph name {" << endl;
	out << "  node [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
	{
		out << "  n" << g.id(n) << " [ label=\"" << index[n] << "\" ]; " << endl;
	}
	out << "  edge [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		out << "  n" << g.id(g.u(e)) << " -- " << " n" << g.id(g.v(e))
				<< " [ label=\"" << length[e] << "\" ]; " << endl;
	}
	out << "}" << endl;

}

typedef pair<unsigned, unsigned> IndexPair;

static void writeEGDot(SmartGraph& g, SmartGraph::EdgeMap<double>& length,
		ostream& out, SmartGraph::NodeMap<IndexPair>& edgePairs,
		unordered_set<unsigned>& matched)
{
	out << "graph name {" << endl;
	out << "  node [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
	{
		IndexPair p = edgePairs[n];
		if(p.first != p.second)
			out << "  n" << g.id(n) << " [ label=\"" << p.first << ":" << p.second << "\" ]; " << endl;
		else
			out << "  n" << g.id(n) << " [ shape=triangle, label=\"" << p.first << "\" ]; \n";
	}
	out << "  edge [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		string extra = "";
		if(matched.count(g.id(e)))
			extra = "style = bold, ";
		out << "  n" << g.id(g.u(e)) << " -- " << " n" << g.id(g.v(e))
				<< " [ " << extra << "label=\"" << length[e] << "\" ]; " << endl;
	}
	out << "}" << endl;
}

/* like knnMergeClusters, but does a quad merge by constructing
 * and edge graph with node constraints from the knn graph
 */
bool MatcherPacker::knnQuadMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache) const
{
	unsigned N = clusters.size();
	SmartGraph SG;
	SmartGraph::EdgeMap<double> Sweights(SG);
	SmartGraph::NodeMap<unsigned> Sindex(SG);
	vector<SmartGraph::Node> Snodes;

	ClusterCache ccache(clusters.size());
	makeKNNGraph(D, clusters, maxSz, dcache, ccache, SG, Sweights, Sindex,
			Snodes);

	SmartGraph::EdgeMap<SmartGraph::Node> edgeNodes(SG);

	SmartGraph EG;
	SmartGraph::NodeMap<IndexPair> edgePairs(EG);
	SmartGraph::EdgeMap<double> Eweights(EG);

	//create node for each edge
	vector<SmartGraph::Node> enodes;
	unsigned enodeCnt = 0;
	for (SmartGraph::EdgeIt e(SG); e != INVALID; ++e)
	{
		SmartGraph::Node en = EG.addNode();
		enodeCnt++;
		edgeNodes[e] = en;
		SmartGraph::Node u = SG.u(e);
		SmartGraph::Node v = SG.v(e);

		edgePairs[en] = IndexPair(Sindex[u], Sindex[v]);
	}

	ofstream sgout("sg.dot");
	writeDot(SG, Sweights, Sindex, sgout);

	unordered_set<IndexPair> seen;

	//now connect edge nodes to edge nodes; the cost is the cost
	//of merging all four shapes
	for (SmartGraph::EdgeIt e(SG); e != INVALID; ++e)
	{
		SmartGraph::Node u = SG.u(e);
		SmartGraph::Node v = SG.v(e);

		//look at all the out edges from u that don't connect back to u or v
		for (SmartGraph::OutArcIt out(SG, u); out != INVALID; ++out)
		{
			SmartGraph::Node u2 = SG.v(out);
			if (u2 != v && u2 != u)
			{
				for (SmartGraph::OutArcIt e2(SG, u2); e2 != INVALID; ++e2)
				{
					SmartGraph::Node v2 = SG.v(e2);
					if (v2 != u && v2 != v && v2 != u2)
					{
						unsigned ei = SG.id(e);
						unsigned e2i = SG.id(e2);

						if (!seen.count(IndexPair(ei, e2i)))
						{
							seen.insert(IndexPair(ei, e2i));
							seen.insert(IndexPair(e2i, ei));

							unsigned a = Sindex[u];
							unsigned b = Sindex[v];
							unsigned c = Sindex[u2];
							unsigned d = Sindex[v2];
							assert(
									a != b && b != c && c != d && a != d && a != c && b != d);

							Eweights[EG.addEdge(edgeNodes[e], edgeNodes[e2])] =
									quadCost(a, b, c, d, clusters, D, ccache,
											dcache);
						}
					}
				}
			}
		}
	}

	seen.clear();

	Mip mip;

	SmartGraph::EdgeMap<Mip::Col> colmap(EG);

	//setup objective function
	Mip::Expr o;
	for (SmartGraph::EdgeIt e(EG); e != INVALID; ++e)
	{
		colmap[e] = mip.addCol();

		cout << "I " << EG.id(e) << " " << mip.id(colmap[e]) << "\n";
		o += Eweights[e]*colmap[e];
		//while we're at it, bound the values
		mip.colLowerBound(colmap[e], 0);
		mip.colUpperBound(colmap[e], 1);

		mip.colType(colmap[e], Mip::INTEGER);
	}
	mip.max();
	mip.obj(o);

	//now constraints, one for each original SG node
	for (SmartGraph::NodeIt n(SG); n != INVALID; ++n)
	{
		unordered_set<unsigned> edgeIds;
		//get all the edges out
		for (SmartGraph::IncEdgeIt sgedge(SG, n); sgedge != INVALID; ++sgedge)
		{
			//each sgedge corresponds to a node in EG, all the outedges
			//of which should be added to the constraint
			SmartGraph::Node en = edgeNodes[sgedge];
			for (SmartGraph::IncEdgeIt egedge(EG, en); egedge != INVALID; ++egedge)
			{
				edgeIds.insert(EG.id(egedge));
			}
		}

		Mip::Expr cons;
		for(unordered_set<unsigned>::iterator itr = edgeIds.begin(), end = edgeIds.end(); itr != end; ++itr)
		{
			cout << *itr << "\t" << mip.id(colmap[EG.edgeFromId(*itr)]) << "\n";
			cons += colmap[EG.edgeFromId(*itr)];
			assert(colmap[EG.edgeFromId(*itr)] != INVALID);
		}

		mip.addRow(cons <= 1);
	}

	mip.solve();

	if(mip.type() == Mip::OPTIMAL)
	{
		cout << "solution " << mip.solValue() << "\n";
	}

	unordered_set<unsigned> setedges;
	for (SmartGraph::EdgeIt e(EG); e != INVALID; ++e)
	{
		if(mip.sol(colmap[e]) > 0)
			setedges.insert(EG.id(e));
	}



	ofstream egout("eg.dot");
	writeEGDot(EG, Eweights, egout, edgePairs, setedges);

	exit(0);

	MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<double> > matcher(EG,
				Eweights);
		matcher.run();
		cout << matcher.matchingWeight() << " " << countNodes(EG) << "\n";


	//merge quads that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 4 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for (unsigned i = 0, n = enodes.size(); i < n; i++)
	{
		SmartGraph::Node e1 = enodes[i];
		unsigned a = edgePairs[e1].first;
		unsigned b = edgePairs[e1].second;

		SmartGraph::Node e2 = matcher.mate(e1);
		if (EG.valid(e2))
		{
			//have a match, all involved better be the same match state
			unsigned c = edgePairs[e2].first;
			unsigned d = edgePairs[e2].second;
			cout << "want to merge " << EG.id(e1) << " " << EG.id(e2) << "\n";
			cout << "want to merge " << a << " " << b << " " << c << " " << d << "\n";
			if(c == d) //a constraint node
				continue;

			if (packed[a] && packed[b] && packed[c] && packed[d])
			{
				//already dealt with
			}
			else if (!packed[a] && !packed[b] && !packed[c] && !packed[d])
			{
				packed[a] = packed[b] = packed[c] = packed[d] = true;
				cout << "merging " << a << " " << b << " " << c << " " << d
						<< "\n";
				if (clusters[a].size() + clusters[b].size() + clusters[c].size()
						+ clusters[d].size() <= packSize)
				{
					//merge
					newclusters.push_back(Cluster());
					newclusters.back().mergeInto(clusters[a], clusters[b],
							clusters[c], clusters[d]);
					didmerge = true;
				}
				else
				{
					//keep clusters as is
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[a]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[b]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[c]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[d]);
				}
			}
			else
			{
				cout << "want to merge " << a << " " << b << " " << c << " "
						<< d << "\n";
				abort();
			}
		}
	}

	//grab anything that wasn't matched
	for (unsigned i = 0; i < N; i++)
	{
		if (!packed[i])
		{
			newclusters.push_back(Cluster());
			newclusters.back().moveInto(clusters[i]);
		}
	}

	swap(clusters, newclusters);
	return didmerge;
}

/* find the optimal matching of every pair of clusters
 * maxSz is what the largest current cluster should be
 */
bool MatcherPacker::fullMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache) const
{
	unsigned C = clusters.size();
	if (C == 1)
		return false;

	//create full graph with distances; enforce an even number of nodes if
	//necessary by adding a dummy node
	unsigned N = C;
	if (C % 2)
		N++;

	FullGraph G(N);
	FullGraph::EdgeMap<double> weights(G);
	for (unsigned i = 0; i < C; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusterDistance(D, clusters[i], clusters[j], dcache);
			weights[G.arc(G(i), G(j))] = -dist; //because the algo maximizes
		}
	}

	//init dummy node
	for (unsigned i = C; i < N; i++)
	{
		for (unsigned j = 0; j < C; j++)
		{
			if (clusters[j].size() < maxSz)
			{
				//force smaller clusters to merge
				weights[G.arc(G(i), G(j))] = -10000;
				weights[G.arc(G(j), G(i))] = -10000;
			}
			else
			{
				weights[G.arc(G(i), G(j))] = 0;
				weights[G.arc(G(j), G(i))] = 0;
			}
		}
	}

	MaxWeightedPerfectMatching<FullGraph, FullGraph::EdgeMap<double> > matcher(
			G, weights);
	matcher.run();

	//merge clusters that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for (unsigned i = 0; i < C; i++)
	{
		if (!packed[i])
		{
			unsigned j = G.index(matcher.mate(G(i)));
			if (j == C)
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
			else
			{
				if (clusters[i].size() + clusters[j].size() <= packSize) //merge
				{
					newclusters.push_back(Cluster());
					newclusters.back().mergeInto(clusters[i], clusters[j]);
					didmerge = true;
					packed[i] = true;
					packed[j] = true;
				}
				else
				{
					//didn't merge, keep clusters as is
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[i]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[j]);
					packed[i] = true;
					packed[j] = true;
				}
			}
		}
	}

	swap(clusters, newclusters);
	return didmerge;
}

//return the cost of merging these four clusters
//two reasonable choices - complete link, where we return the largest distance between any two clusters
//or average linkish - compute the merged MIV/MSV and return the largest distance between any individual and the consensus
double MatcherPacker::quadCost(unsigned a, unsigned b, unsigned c, unsigned d,
		const vector<Cluster>& clusters, const DataViewer *D,
		ClusterCache& cache, DCache& dcache) const
{
	if (true)
	{
		Cluster consensus; //will take control of memory

		const MappableOctTree *itrees[4] =
		{ clusters[a].MIV, clusters[b].MIV, clusters[c].MIV, clusters[d].MIV };

		consensus.MIV = MappableOctTree::createFromIntersection(4, itrees);

		const MappableOctTree *utrees[4] =
		{ clusters[a].MSV, clusters[b].MSV, clusters[c].MSV, clusters[d].MSV };
		consensus.MSV = MappableOctTree::createFromUnion(4, utrees);

		float dist = 0;
		dist = max(dist, clusterDistance(D, clusters[a], consensus, dcache));
		dist = max(dist, clusterDistance(D, clusters[b], consensus, dcache));
		dist = max(dist, clusterDistance(D, clusters[c], consensus, dcache));
		dist = max(dist, clusterDistance(D, clusters[d], consensus, dcache));

		return dist;
	}
	else
	{
		unsigned C[] =
		{ a, b, c, d };
		float dist = 0;
		for (unsigned i = 0; i < 4; i++)
		{
			for (unsigned j = i + 1; j < 4; j++)
			{
				float d = 0;
				if (!cache.isCached(C[i], C[j], d))
				{
					d = clusterDistance(D, clusters[C[i]], clusters[C[j]],
							dcache);
					cache.set(C[i], C[j], d);
				}
				dist = max(dist, d);
			}
		}
		return dist;
	}
}

