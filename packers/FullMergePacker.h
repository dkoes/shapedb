/*
 * FullMergePacker.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *   This packer gready merges all nodes all at once and repeats
 *   until the size limit is hit.
 */

#ifndef FULLMERGEPACKER_H_
#define FULLMERGEPACKER_H_

#include "Packer.h"

class FullMergePacker: public Packer
{
public:

private:
	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters, DCache& dcache) const;

public:

	FullMergePacker(unsigned ps, ClusterDistance metric=AverageLink): Packer(ps, metric) {}
	virtual ~FullMergePacker() {}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* FULLMERGEPACKER_H_ */
