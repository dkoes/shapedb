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
#include <lemon/dimacs.h>

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

//create reverse mapping from forward mapping
static void reverse(const vector<vector<unsigned> >& F,
		vector<vector<unsigned> >& R)
{
	R.clear();
	R.resize(F.size());
	for (unsigned v = 0, n = F.size(); v < n; v++)
	{
		for (unsigned i = 0, m = F[v].size(); i < m; i++)
		{
			unsigned neigh = F[v][i];
			R[neigh].push_back(v);
		}
	}

	for (unsigned v = 0, n = R.size(); v < n; v++)
	{
		sort(R[v].begin(), R[v].end());
	}
}

//insert in sorted order, make K items, no limit if 0
bool MatcherPacker::KNNSlice::update(const IndDist& item, unsigned k)
{
	vector<IndDist>::iterator pos = lower_bound(neighbors.begin(),
			neighbors.end(), item);

	for (unsigned i = 0, n = neighbors.size(); i < n; i++)
	{
		if (neighbors[i].j == item.j && neighbors[i].dist != item.dist)
			abort();
	}
	if (k > 0 && (pos - neighbors.begin()) >= k)
		return false;

	for (vector<IndDist>::iterator look = pos;
			look != neighbors.end() && look->dist == pos->dist; look++)
	{
		if (look->j == item.j)
			return false; //already in here
		if (look->dist == item.dist) //not actually better
			return false;
	}

	neighbors.insert(pos, item);
	if (k > 0)
		neighbors.resize(k);

	return true;
}

//processed and unprocessed; once extracted becomes processed
void MatcherPacker::KNNSlice::getOldNew(vector<unsigned>& O,
		vector<unsigned>& N)
{
	for (unsigned i = 0, n = neighbors.size(); i < n; i++)
	{
		if (neighbors[i].unprocessed)
		{
			N.push_back(neighbors[i].j);
			neighbors[i].unprocessed = false;
		}
		else
			O.push_back(neighbors[i].j);
	}
}

//create an initial knn graph quickly
void MatcherPacker::initialKNNSample(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache,
		ClusterCache& ccache, vector<KNNSlice>& V, float& maxdist) const
{
	if (S == 0) //random sampling
	{
		unsigned N = clusters.size();
		V.resize(N);

		vector<unsigned> indices(N);
		for (unsigned i = 0; i < N; i++)
		{
			indices[i] = i;
		}

		random_shuffle(indices.begin(), indices.end());

		unsigned pos = 0;
		for (unsigned i = 0; i < N; i++)
		{
			V[i].neighbors.reserve(K);
			for (unsigned j = 0; j < K; j++)
			{
				//get random neighbor
				unsigned neigh = indices[pos];
				pos = (pos + 1) % N;
				if (neigh == i)
				{
					neigh = indices[pos];
					pos = (pos + 1) % N;
				}

				float dist = 0;
				if (!ccache.isCached(i, neigh, dist))
				{
					dist = clusterDistance(D, clusters[i], clusters[neigh],
							dcache);
					ccache.set(i, neigh, dist);
				}

				if (dist > maxdist)
					maxdist = dist;
				V[i].neighbors.push_back(IndDist(neigh, dist));
			}

			sort(V[i].neighbors.begin(), V[i].neighbors.end());
		}
	}
	else
	{
		//pick S sentinals
		unsigned N = clusters.size();
		V.resize(N);

		ANNpointArray points = annAllocPts(N, S);

		//farthest first, initialize point data as we go
		vector<unsigned> indices;
		indices.push_back(0);
		for (unsigned s = 1; s < S; s++)
		{
			float max = 0;
			unsigned besti = 0;
			for (unsigned i = 0; i < N; i++)
			{
				float min = HUGE_VAL;
				unsigned send = indices.size() - 1; //all but last have distances already computed
				for (unsigned j = 0; j < send; j++)
				{
					float dist = points[i][j];
					if (dist < min)
						min = dist;
				}
				float dist = 0;
				if (i != indices.back())
				{
					if (!ccache.isCached(i, indices.back(), dist))
					{
						dist = clusterDistance(D, clusters[i],
								clusters[indices.back()], dcache);
						ccache.set(i, indices.back(), dist);
					}
				}
				points[i][send] = dist;
				if (dist < min)
					min = dist;

				if (min > max)
				{
					besti = i;
					max = min;
				}
			}
			indices.push_back(besti);
		}assert(indices.size() == S);
		//get distances for last chosen sentinal
		for (unsigned i = 0; i < N; i++)
		{
			float dist = 0;
			if (i != indices.back())
			{
				if (!ccache.isCached(i, indices.back(), dist))
				{
					dist = clusterDistance(D, clusters[i],
							clusters[indices.back()], dcache);
					ccache.set(i, indices.back(), dist);
				}
			}
			points[i][S - 1] = dist;
		}

		ANNkd_tree searcher(points, N, S);

		ANNdistArray dists = new ANNdist[K + 1];
		vector<int> nnIdx(K + 1, 0);
		for (unsigned i = 0; i < N; i++)
		{
			V[i].neighbors.reserve(K);
			searcher.annkSearch(points[i], K + 1, &nnIdx[0], dists);

			for (unsigned j = 0; j <= K; j++)
			{
				if (nnIdx[j] != (int) i && nnIdx[j] >= 0)
				{
					unsigned neigh = nnIdx[j];
					float dist = 0;
					if (!ccache.isCached(i, neigh, dist))
					{
						dist = clusterDistance(D, clusters[i], clusters[neigh],
								dcache);
						ccache.set(i, neigh, dist);
					}

					V[i].neighbors.push_back(IndDist(nnIdx[j], dist));

					if (dist > maxdist)
						maxdist = dist;
				}
			}
		}

		delete[] dists;
		annDeallocPts(points);
	}
}

//create a knn graph, return max value of any edge since we reverse the meanings of distances
float MatcherPacker::makeKNNGraph(const DataViewer *D, vector<Cluster>& clusters,
		unsigned maxsz, DCache& dcache, ClusterCache& ccache, SmartGraph& G,
		SmartGraph::EdgeMap<double>& E, SmartGraph::NodeMap<unsigned>& Sindex,
		vector<SmartGraph::Node>& nodes) const
{
	G.clear();
	unsigned N = clusters.size();
	float max = 0;

	if (K >= clusters.size()) //full knn
	{
		vector<vector<IndDist> > dists(N, vector<IndDist>(N));
		for (unsigned i = 0; i < N; i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				float dist = 0;

				if (!ccache.isCached(i, j, dist))
				{
					dist = clusterDistance(D, clusters[i], clusters[j], dcache);
					ccache.set(i, j, dist);
				}
				unsigned minsz = min(clusters[i].size(), clusters[j].size());
				while (minsz < maxsz)
				{
					dist /= 2;
					minsz *= 2;
				}
				dists[i][j] = IndDist(j, dist);
				dists[j][i] = IndDist(i, dist);

				if (dist > max)
					max = dist;
			}
			dists[i][i] = IndDist(i, HUGE_VAL);
		}

		G.reserveNode(N);
		G.reserveEdge(N * K);

		nodes.resize(N);
		for (unsigned i = 0; i < N; i++)
		{
			nodes[i] = G.addNode();
			Sindex[nodes[i]] = i;
		}

		for (unsigned i = 0; i < N; i++)
		{
			sort(dists[i].begin(), dists[i].end());

			for (unsigned j = 0; j < K && j < N; j++)
			{
				if (dists[i][j].j != i)
				{
					SmartGraph::Edge e = G.addEdge(nodes[i],
							nodes[dists[i][j].j]);
					E[e] = max - dists[i][j].dist; //compute max weighting
				}
			}
		}
	}
	else //approximate knn
	{
		vector<KNNSlice> B(N);

		initialKNNSample(D, clusters, maxsz, dcache, ccache, B, max);

		bool keepgoing = true;
		while (keepgoing)
		{
			keepgoing = false;
			unsigned changed = 0;
			vector<vector<unsigned> > oldf(N), newf(N);
			vector<vector<unsigned> > oldr(N), newr(N);
			for (unsigned v = 0; v < N; v++)
			{
				B[v].getOldNew(oldf[v], newf[v]);
			}

			reverse(oldf, oldr);
			reverse(newf, newr);

			for (unsigned v = 0; v < N; v++)
			{
				//merge forward and reverse neighbors
				vector<unsigned> oldn, newn;
				sort(oldf[v].begin(), oldf[v].end());
				sort(oldr[v].begin(), oldr[v].end());
				sort(newf[v].begin(), newf[v].end());
				sort(newr[v].begin(), newr[v].end());

				merge(oldf[v].begin(), oldf[v].end(), oldr[v].begin(),
						oldr[v].end(), back_inserter(oldn));
				merge(newf[v].begin(), newf[v].end(), newr[v].begin(),
						newr[v].end(), back_inserter(newn));

				//compare new to new and new to old
				for (unsigned i = 0, n = newn.size(); i < n; i++)
				{
					unsigned u1 = newn[i];
					for (unsigned j = 0; j < n; j++)
					{
						unsigned u2 = newn[j];
						if (u1 < u2)
						{
							float dist = 0;
							if (!ccache.isCached(u1, u2, dist))
							{
								dist = clusterDistance(D, clusters[u1],
										clusters[u2], dcache);
								ccache.set(u1, u2, dist);
							}

							changed += B[u1].update(IndDist(u2, dist), K);
							changed += B[u2].update(IndDist(u1, dist), K);

							if (dist > max)
								max = dist;
						}
					}

					for (unsigned j = 0, m = oldn.size(); j < m; j++)
					{
						unsigned u2 = oldn[j];
						if (u2 != u1)
						{
							float dist = 0;
							if (!ccache.isCached(u1, u2, dist))
							{
								dist = clusterDistance(D, clusters[u1],
										clusters[u2], dcache);
								ccache.set(u1, u2, dist);
							}

							if (B[u1].update(IndDist(u2, dist), K))
								changed++;
							if (B[u2].update(IndDist(u1, dist), K))
								changed++;

							if (dist > max)
								max = dist;
						}
					}
				}
			}

			if (changed > 0)
				keepgoing = true;
		}

		G.reserveNode(N);
		G.reserveEdge(N * K);

		nodes.resize(N);
		for (unsigned i = 0; i < N; i++)
		{
			nodes[i] = G.addNode();
			Sindex[nodes[i]] = i;
		}

		typedef pair<unsigned, unsigned> UPair;
		unordered_set<UPair> seen;
		for (unsigned i = 0; i < N; i++)
		{
			for (unsigned j = 0, nn = B[i].neighbors.size(); j < nn; j++)
			{
				unsigned neigh = B[i].neighbors[j].j;

				if (seen.count(UPair(i, neigh)) == 0)
				{
					seen.insert(UPair(i, neigh));
					seen.insert(UPair(neigh, i));
					float dist = B[i].neighbors[j].dist;
					SmartGraph::Edge e = G.addEdge(nodes[i], nodes[neigh]);
					E[e] = max - dist; //compute max weighting
				}
			}
		}
	}
	return max;
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
		MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<double> >& matcher)
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
		if(matcher.matching(e))
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
	float max = makeKNNGraph(D, clusters, maxSz, dcache, ccache, SG, Sweights, Sindex,
			Snodes);

	OutDegMap<SmartGraph> degrees(SG);
	SmartGraph::NodeMap<vector<SmartGraph::Node> > constrainNodes(SG);
	SmartGraph::EdgeMap<SmartGraph::Node> edgeNodes(SG);

	SmartGraph EG;
	SmartGraph::NodeMap<IndexPair> edgePairs(EG);
	SmartGraph::EdgeMap<double> Eweights(EG);
	//construct the constrained edge graph - every edge becomes a node, plus
	//node constrained required degreee-1 additional nodes for each node
	unsigned cnodeCnt = 0;
	vector<SmartGraph::Node> allcnodes;
	for (SmartGraph::NodeIt n(SG); n != INVALID; ++n)
	{
		unsigned deg = degrees[n];
		if (deg > 0)
			deg--;
		vector<SmartGraph::Node> cnodes;
		cnodes.reserve(deg);
		for (unsigned i = 0; i < deg; i++)
		{
			cnodes.push_back(EG.addNode());
			allcnodes.push_back(cnodes.back());
			cnodeCnt++;
			edgePairs[cnodes.back()] = IndexPair(Sindex[n], Sindex[n]); //for debug, rm
		}
		constrainNodes[n] = cnodes;
	}


	//create node for each edge
	vector<SmartGraph::Node> enodes;
	unsigned enodeCnt = 0;
	for (SmartGraph::EdgeIt e(SG); e != INVALID; ++e)
	{
		SmartGraph::Node en = EG.addNode();
		enodeCnt++;
		edgeNodes[e] = en;
		//conect to constraintNodes
		SmartGraph::Node u = SG.u(e);
		SmartGraph::Node v = SG.v(e);

		enodes.push_back(en);
		edgePairs[en] = IndexPair(Sindex[u], Sindex[v]);

		const vector<SmartGraph::Node> unodes = constrainNodes[u];
		for (unsigned i = 0, n = unodes.size(); i < n; i++)
		{
			Eweights[EG.addEdge(unodes[i], en)] = 2*max;
		}
		const vector<SmartGraph::Node> vnodes = constrainNodes[v];
		for (unsigned i = 0, n = vnodes.size(); i < n; i++)
		{
			Eweights[EG.addEdge(vnodes[i], en)] = 2*max;
		}
	}

	ofstream sgout("sg.dot");
	writeDot(SG, Sweights, Sindex, sgout);

	unordered_set<IndexPair> seen;

	//now connect edge nodes to edge nodes; the cost is the cost
	//of merging all four shapes
	for (SmartGraph::EdgeIt e(SG); e != INVALID; ++e)
	{
		//conect to constraintNodes
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

	MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<double> > matcher(EG,
			Eweights);
	matcher.run();
	cout << matcher.matchingWeight() << " " << countNodes(EG) << "\n";

	ofstream egout("eg.dot");
	writeEGDot(EG, Eweights, egout, edgePairs, matcher);


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

