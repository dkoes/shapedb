/*
 * Packer.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef PACKER_H_
#define PACKER_H_

#include "DataViewers.h"
#include "ShapeDistance.h"
#include <boost/multi_array.hpp>


//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
public:
	enum ClusterDistance { CompleteLink, AverageLink, SingleLink, TotalLink, NotApplicable };
protected:
	unsigned packSize;
	ClusterDistance distMetric;

	//caches distance information between single shapes
	class DCache
	{
		multi_array<float, 2> cache;

	public:
		//indices are assumed to be consecutive and sorted
		DCache(const DataViewer* D): cache(extents[D->size()][D->size()])
		{
			unsigned N = D->size();
			for(unsigned i = 0; i < N; i++)
			{
				for(unsigned j = 0; j < i; j++)
				{
					float dist = shapeDistance(D->getMIV(i), D->getMSV(i), D->getMIV(j),
							D->getMSV(j));
					cache[i][j] = dist;
					cache[j][i] = dist;
				}
				cache[i][i] = 0;
			}
		}

		float get(unsigned i, unsigned j)
		{
			return cache[i][j];
		}
	};

	float clusterDistance(const DataViewer* D, const Cluster& a,
			const Cluster& b, DCache& dcache) const;
public:

	Packer(unsigned ps, ClusterDistance metric): packSize(ps), distMetric(metric) {}
	virtual ~Packer() {}
	//assumes dv has been sliced into sequential indices
	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const = 0;

	unsigned getPack() const { return packSize; }
};



#endif /* PACKER_H_ */
