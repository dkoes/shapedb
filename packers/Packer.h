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
	enum ClusterDistance
	{
		CompleteLink, AverageLink, SingleLink, TotalLink, NotApplicable
	};
protected:
	unsigned packSize;
	ClusterDistance distMetric;

	//caches distance information between single shapes
	class DCache
	{
		multi_array<float, 2> cache;
		const DataViewer* dv;

	public:
		//indices are assumed to be consecutive and sorted
		DCache(const DataViewer* D, bool preload = true) :
				cache(extents[D->size()][D->size()]), dv(D)
		{
			if (preload)
			{
				unsigned N = D->size();
				for (unsigned i = 0; i < N; i++)
				{
					for (unsigned j = 0; j < i; j++)
					{
						float dist = shapeDistance(D->getMIV(i), D->getMSV(i),
								D->getMIV(j), D->getMSV(j));
						cache[i][j] = dist;
						cache[j][i] = dist;
					}
					cache[i][i] = 0;
				}
			}
			else
			{
				unsigned N = D->size();
				for (unsigned i = 0; i < N; i++)
				{
					for (unsigned j = 0; j < N; j++)
					{
						cache[i][j] = -HUGE_VAL;
					}
				}
			}
		}

		float get(unsigned i, unsigned j)
		{
			if (cache[i][j] < 0)
			{
				if (i == j)
					cache[i][j] = 0;
				else
				{
					float dist = shapeDistance(dv->getMIV(i), dv->getMSV(i),
							dv->getMIV(j), dv->getMSV(j));
					cache[i][j] = dist;
					cache[j][i] = dist;
				}
			}
			return cache[i][j];
		}

		unsigned cnt() const
		{
			unsigned N = cache.size();
			unsigned ret = 0;
			for (unsigned i = 0; i < N; i++)
			{
				for (unsigned j = 0; j < i; j++)
				{
					if(cache[i][j] >= 0)
						ret++;
				}
			}
			return ret;
		}
	};

	float clusterDistance(const DataViewer* D, const Cluster& a,
			const Cluster& b, DCache& dcache) const;
public:

	Packer(unsigned ps, ClusterDistance metric) :
			packSize(ps), distMetric(metric)
	{
	}
	virtual ~Packer()
	{
	}
	//assumes dv has been sliced into sequential indices
	virtual void pack(const DataViewer* dv,
			vector<Cluster>& clusters) const = 0;

	unsigned getPack() const
	{
		return packSize;
	}
};

#endif /* PACKER_H_ */
