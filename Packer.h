/*
 * Packer.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef PACKER_H_
#define PACKER_H_

#include "DataViewers.h"
//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
public:
	enum ClusterDistance { CompleteLink, AverageLink, SingleLink, NotApplicable };
protected:
	unsigned packSize;
	ClusterDistance distMetric;
public:

	Packer(unsigned ps, ClusterDistance metric): packSize(ps), distMetric(metric) {}
	virtual ~Packer() {}
	virtual void pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const = 0;

	unsigned getPack() const { return packSize; }
};



#endif /* PACKER_H_ */
