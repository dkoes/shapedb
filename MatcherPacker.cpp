/*
 * MatcherPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "MatcherPacker.h"
#include <lemon/full_graph.h>
#include <lemon/matching.h>

using namespace lemon;

void MatcherPacker::pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const
{
	if (indices.size() == 0)
		return;

	//start with each object in its own cluster
	clusters.clear();
	clusters.resize(indices.size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	//combine everything as much as possible
	while (fullMergeClusters(dv, clusters))
		;
}


bool MatcherPacker::fullMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters) const
{
	unsigned C = clusters.size();
	if (C == 1)
		return false;

	//create full graph with distances; enforce an even number of nodes if
	//necessary by adding a dummy node
	unsigned N = C;
	if(C%2)
		N++;
	FullGraph G(N);
	FullGraph::EdgeMap<double> weights(G);
	for (unsigned i = 0; i < C; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusterDistance(D, clusters[i], clusters[j]);
			if(clusters[i].size() + clusters[j].size() > packSize)
				dist = 1000000;
			weights[G.arc(G(i),G(j))] = -dist; //because the algo maximizes
		}
	}

	//init dummy node
	for(unsigned i = C; i < N; i++)
	{
		for(unsigned j = 0; j < C; j++)
		{
			weights[G.arc(G(i),G(j))] = 0;
			weights[G.arc(G(j),G(i))] = 0;
		}
	}

	MaxWeightedPerfectMatching<FullGraph, FullGraph::EdgeMap<double>  > matcher(G,weights);
	matcher.run();

	//merge clusters that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for(unsigned i = 0; i < C; i++)
	{
		if(!packed[i])
		{
			unsigned j = G.index(matcher.mate(G(i)));
			if(j == C) //dummy
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
			else if(clusters[i].size() + clusters[j].size() <= packSize)//merge
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

	swap(clusters,newclusters);
	return didmerge;
}
