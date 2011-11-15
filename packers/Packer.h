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

	class DCache //interface for caching distance information between shapes
	{
	public:
		virtual float get(unsigned i, unsigned j) = 0;

	};

	//stores and precomputes everything
	class FullCache: public DCache
	{
		multi_array<float, 2> cache;
		const DataViewer* dv;

	public:
		//indices are assumed to be consecutive and sorted
		FullCache(const DataViewer* D) :
				cache(extents[D->size()][D->size()]), dv(D)
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

		float get(unsigned i, unsigned j)
		{
			return cache[i][j];
		}

		unsigned size() const
		{
			return cache.size() * cache.size();
		}
	};

	//does not precompute, uses much less memory
	class OnDemandCache: public DCache
	{
		struct JDist
		{
			unsigned j;
			float dist;

			JDist(): j(0), dist(0) {}
			JDist(unsigned J, float D): j(J), dist(D) {}

			bool operator<(const JDist& rhs) const
			{
				return j < rhs.j;
			}

			bool operator==(const JDist& rhs) const
			{
				return j == rhs.j;
			}
		};

		vector< vector<JDist> > cache;
		const DataViewer* dv;

	public:
		//indices are assumed to be consecutive and sorted
		OnDemandCache(const DataViewer* D) :
				cache(D->size()), dv(D)
		{

		}

		float get(unsigned i, unsigned j)
		{
			if(j > i) //normalize order - store only once
				swap(i,j);
			else if(i == j)
				return 0;

			vector<JDist>& row = cache[i];
			JDist val(j,0);
			vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(), val);

			if(pos != row.end() && pos->j == j)
			{
				//already exists
				return pos->dist;
			}
			else
			{
				//must compute and insert
				val.dist = shapeDistance(dv->getMIV(i), dv->getMSV(i),
						dv->getMIV(j), dv->getMSV(j));
				row.insert(pos, val);
				return val.dist;
			}
		}

		unsigned size() const
		{
			unsigned ret = 0;
			for(unsigned i = 0, n = cache.size(); i < n; i++)
			{
				ret += cache[i].size();
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
