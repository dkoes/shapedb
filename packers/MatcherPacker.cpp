/*
 * MatcherPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "MatcherPacker.h"
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>

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
		DCache dcache(dv);

		while (curSz < packSize
				&& fullMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;
	}
	else //use knn
	{
		DCache dcache(dv, false);

		while (curSz < packSize && knnMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;
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

	if (k > 0 && (pos - neighbors.begin()) >= k)
		return false;

	for(vector<IndDist>::iterator look = pos; look != neighbors.end() && look->dist == pos->dist; look++)
	{
		if(look->j == item.j)
			return false; //already in here
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
		if (neighbors[i].unprocessed) {
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
		vector<KNNSlice>& V) const
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

			float dist = clusterDistance(D, clusters[i], clusters[neigh],
					dcache);
			V[i].neighbors.push_back(IndDist(neigh, dist));
		}

		sort(V[i].neighbors.begin(), V[i].neighbors.end());
	}
}

//create a knn graph
void MatcherPacker::makeKNNGraph(const DataViewer *D, vector<Cluster>& clusters,
		unsigned maxsz, DCache& dcache, SmartGraph& G,
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
				float dist = clusterDistance(D, clusters[i], clusters[j],
						dcache);
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

		initialKNNSample(D, clusters, maxsz, dcache, B);

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
							float dist = clusterDistance(D, clusters[u1],
									clusters[u2], dcache);
							changed += B[u1].update(IndDist(u2, dist), K);
							changed += B[u2].update(IndDist(u1, dist), K);
						}
					}

					for (unsigned j = 0, m = oldn.size(); j < m; j++)
					{
						unsigned u2 = oldn[j];
						float dist = clusterDistance(D, clusters[u1],
								clusters[u2], dcache);
						changed += B[u1].update(IndDist(u2, dist), K);
						changed += B[u2].update(IndDist(u1, dist), K);
					}
				}
			}

			if (changed > 0)
				keepgoing = true;
			cout << changed << "\n";
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
			for (unsigned j = 0, nn = B[i].neighbors.size(); j < nn; j++)
			{
				unsigned neigh = B[i].neighbors[j].j;
				float dist = B[i].neighbors[j].dist;
				SmartGraph::Edge e = G.addEdge(nodes[i], nodes[neigh]);
				E[e] = 10 - dist; //compute max weighting
			}
		}
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
	makeKNNGraph(D, clusters, maxSz, dcache, SG, Sweights, Sindex, Snodes);

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
