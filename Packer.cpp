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
		const Cluster& b) const
{
	switch (distMetric)
	{
	case AverageLink:
		return shapeDistance(a.MIV, a.MSV, b.MIV, b.MSV);
	case CompleteLink:
	{
		//this is the maximum of the minimum distances between cluster members
		//TODO: this is horribly inefficient due to distance recomputation
		float max = 0;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			float min = HUGE_VAL;
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				float dist = 0;
				unsigned l = a[i];
				unsigned r = b[j];

				dist = shapeDistance(D->getMIV(l), D->getMSV(l), D->getMIV(r),
						D->getMSV(r));

				if (dist < min)
					min = dist;
			}
			if (min > max)
				max = min;
		}
		return max;
	}
	case SingleLink:
	{
		//the minimum distance overall
		float min = HUGE_VAL;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				float dist = 0;

				unsigned l = a[i];
				unsigned r = b[j];
				dist = shapeDistance(D->getMIV(l), D->getMSV(l), D->getMIV(r),
						D->getMSV(r));

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
