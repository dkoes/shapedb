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

		void getCommonNeighbors(const vector<unsigned>& rev, vector<unsigned>& neighs) const;
	};

	void initialKNNSample(const DataViewer *D,
			vector<Cluster>& clusters, unsigned maxSz, DCache& dcache, vector< KNNSlice >& V) const;

	void makeKNNGraph(const DataViewer *D,
			vector<Cluster>& clusters, unsigned maxSz, DCache& dcache, lemon::SmartGraph& G,
			lemon::SmartGraph::EdgeMap<double>& E,
			lemon::SmartGraph::NodeMap<unsigned>& Sindex,
			vector<lemon::SmartGraph::Node>& nodes) const;

	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	bool knnMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	unsigned K; //for knn, if zero use full
	unsigned S; //number of sentinals for knn building, if zero, do random
public:
	MatcherPacker(unsigned ps, unsigned k=0, unsigned s=0, ClusterDistance metric = AverageLink) :
			Packer(ps, metric), K(k), S(s)
	{
	}
	~MatcherPacker()
	{
	}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* MATCHERPACKER_H_ */
