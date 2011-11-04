/*
 * GreedyPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "GreedyPacker.h"
#include <boost/multi_array.hpp>

//distances between i and j
struct IntraClusterDist
{
	unsigned i;
	unsigned j;
	float dist;

	IntraClusterDist(): i(0), j(0), dist(HUGE_VAL) {}

	IntraClusterDist(unsigned I, unsigned J, float d): i(I), j(J), dist(d) {}

	bool operator<(const IntraClusterDist& rhs) const
	{
		return dist < rhs.dist;
	}
};


typedef bool (*ICDComp)(const IntraClusterDist& lhs, const IntraClusterDist& rhs);
static bool reverseDist(const IntraClusterDist& lhs, const IntraClusterDist& rhs)
{
	return lhs.dist > rhs.dist;
}

//a priority q for distances, roll our own to support more efficient
//removal of clusters and re-use of distances array
class DistancePQ
{
	vector<IntraClusterDist>& Q;

public:
	DistancePQ(vector<IntraClusterDist>& q): Q(q)
	{
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

	bool empty() const { return Q.size() == 0; }

	const IntraClusterDist& top() const { return Q[0]; }

	void push(const IntraClusterDist& c)
	{
		Q.push_back(c);
		push_heap(Q.begin(), Q.end(), reverseDist);
	}

	void pop()
	{
		pop_heap(Q.begin(), Q.end(), reverseDist);
		Q.pop_back();
	}

	//remove all references to clusters i and j; this is linear
	void removeIJ(unsigned i, unsigned j)
	{
		unsigned c = 0;
		unsigned end = Q.size();
		while(c < end)
		{
			if(Q[c].i != i && Q[c].i != j && Q[c].j != i && Q[c].j != j)
			{
				c++;
			}
			else
			{
				end--;
				swap(Q[c],Q[end]);
			}
		}

		Q.erase(Q.begin()+end, Q.end());
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

};

void GreedyPacker::pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& finalclusters) const
{
	if (indices.size() == 0)
		return;

	//start with each object in its own cluster
	unsigned N = indices.size();
	vector<Cluster> clusters;
	clusters.reserve(N*2);
	clusters.resize(N);
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	//compute all intra cluster distances
	vector<IntraClusterDist> distances; distances.reserve(N*N / 2);
	boost::multi_array<float, 2> darray(boost::extents[N][N]);
	DCache dcache;

	for (unsigned i = 0; i < N; i++)
	{
		darray[i][i] = 0;
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusterDistance(dv, clusters[i], clusters[j], dcache);
			darray[i][j] = darray[j][i] = dist;
			distances.push_back(IntraClusterDist(i,j,dist));
		}
	}

	DistancePQ pQ(distances); //manipulates a reference of distances

	unsigned max = packSize;
	while(!pQ.empty())
	{
		//get the smallest distance
		unsigned i = pQ.top().i;
		unsigned j = pQ.top().j;
		pQ.pop();

		if(clusters[i].isValid() && clusters[j].isValid())
		{
			if(clusters[i].size() + clusters[j].size() <= max)
			{
				unsigned c = clusters.size();
				//merge, create a new cluster
				clusters.push_back(Cluster());
				clusters.back().mergeInto(clusters[i], clusters[j]);

				pQ.removeIJ(i, j);
				//now compute the distance between this cluster and all remaining clusters
				//and add to priority Q
				for(unsigned d = 0; d < c; d++)
				{
					if(clusters[d].isValid())
					{
						float dist = clusterDistance(dv, clusters.back(), clusters[d], dcache);
						pQ.push(IntraClusterDist(c,d,dist));
					}
				}
			}
		}
	}

	//remove empty clusters
	finalclusters.clear();
	finalclusters.reserve(ceil(N/(double)packSize));
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		if(clusters[i].isValid())
		{
			finalclusters.push_back(Cluster());
			finalclusters.back().moveInto(clusters[i]);
		}
	}
}
