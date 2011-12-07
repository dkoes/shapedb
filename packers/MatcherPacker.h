/*
 * MatcherPacker.h
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 *
 *      This packer is similar to the FullMerge packer, but computes the optimal
 *      edge matching at each merge.
 */

#ifndef MATCHERPACKER_H_
#define MATCHERPACKER_H_

#include "Packer.h"
#include <lemon/smart_graph.h>

//caches intercluster distances
class ClusterCache
{
	struct JDist
	{
		unsigned j;
		float dist;

		JDist() :
				j(0), dist(0)
		{
		}
		JDist(unsigned J, float D) :
				j(J), dist(D)
		{
		}

		bool operator<(const JDist& rhs) const
		{
			return j < rhs.j;
		}

		bool operator==(const JDist& rhs) const
		{
			return j == rhs.j;
		}
	};

	vector<vector<JDist> > cache;

public:
	//indices are assumed to be consecutive and sorted
	ClusterCache(unsigned n) :
			cache(n)
	{

	}

	bool isCached(unsigned i, unsigned j, float& ret)
	{
		if (j > i) //normalize order - store only once
			swap(i, j);
		else if (i == j)
		{
			ret = 0;
			return true;
		}

		if (i >= cache.size())
			return false;

		vector<JDist>& row = cache[i];
		JDist val(j, 0);
		vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(), val);

		if (pos != row.end() && pos->j == j)
		{
			//already exists
			ret = pos->dist;
			return true;
		}
		else
		{
			return false;
		}
	}

	void set(unsigned i, unsigned j, float v)
	{
		if (j > i) //normalize order - store only once
			swap(i, j);
		else if (i == j)
			return;

		if (i >= cache.size())
			cache.resize(i + 1);

		vector<JDist>& row = cache[i];
		JDist val(j, v);
		vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(), val);

		if (pos != row.end() && pos->j == j)
		{
			//already exists
			pos->dist = v;
		}
		else
		{
			row.insert(pos, val);
		}
	}
};

class MatcherPacker: public Packer
{
	struct IndDist
	{
		unsigned j;
		bool unprocessed;
		double dist;

		IndDist() :
				j(0), unprocessed(true), dist(HUGE_VAL)
		{
		}
		IndDist(unsigned J, double D) :
				j(J), unprocessed(true), dist(D)
		{
		}

		bool operator<(const IndDist& rhs) const
		{
			return dist < rhs.dist;
		}
	};

	struct KNNSlice
	{
		vector<IndDist> neighbors;

		bool update(const IndDist& item, unsigned k); //insert in sorted order, make K items, no limit if 0
		void getOldNew(vector<unsigned>& o, vector<unsigned>& n); //processed and unprocessed; once extracted becomes processed

		void getCommonNeighbors(const vector<unsigned>& rev,
				vector<unsigned>& neighs) const;
	};

	void initialKNNSample(const DataViewer *D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache, ClusterCache& ccache,
			vector<KNNSlice>& V, float& maxdist) const;

	void makeKNNGraph(const DataViewer *D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache, ClusterCache& ccache,
			lemon::SmartGraph& G, lemon::SmartGraph::EdgeMap<double>& E,
			lemon::SmartGraph::NodeMap<unsigned>& Sindex,
			vector<lemon::SmartGraph::Node>& nodes) const;

	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	bool knnMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	double quadCost(unsigned a, unsigned b, unsigned c, unsigned d,
			const vector<Cluster>& clusters, const DataViewer* D,
			ClusterCache& cache, DCache& dcache) const;

	bool knnQuadMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	unsigned K; //for knn, if zero use full
	unsigned S; //number of sentinals for knn building, if zero, do random
	bool doQuadPack; //if true combine quads
public:
	MatcherPacker(unsigned ps, unsigned k = 0, unsigned s = 0,
			ClusterDistance metric = AverageLink, bool doQ = false) :
			Packer(ps, metric), K(k), S(s), doQuadPack(doQ)
	{
	}
	~MatcherPacker()
	{
	}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* MATCHERPACKER_H_ */
