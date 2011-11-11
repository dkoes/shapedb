/*
 * Packer.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "Packer.h"
#include "ShapeDistance.h"


//return the distance between two clusters, this may be configurable to other metrics
float Packer::clusterDistance(const DataViewer* D, const Cluster& a,
		const Cluster& b, DCache& dcache) const
{
	switch (distMetric)
	{
	case AverageLink:
		if(a.size() == 1 && b.size() == 1)
		{
			return dcache.get(a[0],b[0]);
		}
		else
		{
			return shapeDistance(a.MIV, a.MSV, b.MIV, b.MSV);
		}
	case CompleteLink:
	{
		//this is the maximum o distance between any two cluster members
		float max = 0;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				if (dist > max)
					max = dist;
			}
		}
		return max;
	}
	case TotalLink:
	{
		//this is the total sum of all the linkages, distiguished from average link
		//in that smaller clusters end up with smaller values
		float sum = 0;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				sum += dist;
			}
		}
		return sum;
	}
	case SingleLink:
	{
		//the minimum distance overall
		float min = HUGE_VAL;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				if (dist < min)
					min = dist;
			}
		}
		return min;
	}
	default:
		abort();
		break;
	}
	return 0;
}
