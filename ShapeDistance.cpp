/*
 * ShapeDistance.cpp
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 *
 *      Defines various metrics and initializes shapeDistance to one of them
 */


#include "ShapeDistance.h"


float averageRelVolume(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->relativeVolumeDistance(rightMIV);
	}

	return leftMIV->relativeVolumeDistance(rightMIV) + leftMSV->relativeVolumeDistance(rightMSV);
}

float averageAbsVolume(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->absoluteVolumeDistance(rightMIV);
	}

	return leftMIV->absoluteVolumeDistance(rightMIV) + leftMSV->relativeVolumeDistance(rightMSV);
}

float averageRelVolumeSim(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return 1-leftMIV->relativeVolumeDistance(rightMIV);
	}

	return ((1-leftMIV->relativeVolumeDistance(rightMIV)) + (1-leftMSV->relativeVolumeDistance(rightMSV)))/2.0;
}

shapeMetricFn shapeDistance = averageRelVolume;
shapeMetricFn shapeSimilarity = averageRelVolumeSim;

