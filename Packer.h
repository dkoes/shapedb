/*
 * Packer.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef PACKER_H_
#define PACKER_H_

#include "DataViewers.h"
#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>

//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
public:
	typedef pair<unsigned,unsigned> DCacheKey;
	typedef unordered_map< DCacheKey, double> DCache;
	enum ClusterDistance { CompleteLink, AverageLink, SingleLink, NotApplicable };
protected:
	unsigned packSize;
	ClusterDistance distMetric;

	float clusterDistance(const DataViewer* D, const Cluster& a,
			const Cluster& b, DCache& dcache) const;
public:

	Packer(unsigned ps, ClusterDistance metric): packSize(ps), distMetric(metric) {}
	virtual ~Packer() {}
	virtual void pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const = 0;

	unsigned getPack() const { return packSize; }
};



#endif /* PACKER_H_ */
