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


class MatcherPacker: public Packer
{
	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	bool knnMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	double quadCost(unsigned a, unsigned b, unsigned c, unsigned d,
			const vector<Cluster>& clusters, const DataViewer* D,
			ClusterCache& cache, DCache& dcache) const;

	bool knnQuadMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;


	bool doQuadPack; //if true combine quads EXPERIMENTAL
public:
	MatcherPacker(unsigned ps, unsigned k = 0, unsigned s = 0,
			ClusterDistance metric = AverageLink, bool doQ = false) :
			Packer(ps, metric, k, s), doQuadPack(doQ)
	{
	}
	~MatcherPacker()
	{
	}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* MATCHERPACKER_H_ */
