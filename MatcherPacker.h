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
	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters, DCache& dcache) const;

public:
	MatcherPacker(unsigned ps, ClusterDistance metric=AverageLink): Packer(ps, metric) {}
	~MatcherPacker() {}

	virtual void pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const;

};

#endif /* MATCHERPACKER_H_ */
