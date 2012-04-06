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

//this combines MIV/MSV similarity with the tightness of the resulting MIV/MSV
float tripleRelative(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->relativeVolumeDistance(rightMIV); //highly preferential to merging leaves
	}

	float mivIntersect = 0, mivUnion = 0, msvIntersect = 0, msvUnion = 0;
	leftMIV->intersectUnionVolume(rightMIV, mivIntersect, mivUnion);
	leftMSV->intersectUnionVolume(rightMSV, msvIntersect, msvUnion);

	float mivVol =  1 - mivIntersect/mivUnion;
	float msvVol = 1 - msvIntersect/msvUnion;
	float resultVol = 1 - mivIntersect/msvUnion;
	return mivVol+msvVol+resultVol;
}

//this combines MIV/MSV similarity with the tightness of the resulting MIV/MSV
float tripleAbsolute(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->absoluteVolumeDistance(rightMIV);
	}

	float mivIntersect = 0, mivUnion = 0, msvIntersect = 0, msvUnion = 0;
	leftMIV->intersectUnionVolume(rightMIV, mivIntersect, mivUnion);
	leftMSV->intersectUnionVolume(rightMSV, msvIntersect, msvUnion);

	float mivVol =  mivUnion-mivIntersect;
	float msvVol = msvUnion-msvIntersect;
	float resultVol = msvUnion-mivIntersect;
	return (mivVol+msvVol+resultVol)/3.0;
}

float averageAbsVolume(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->absoluteVolumeDistance(rightMIV);
	}

	return leftMIV->absoluteVolumeDistance(rightMIV) + leftMSV->absoluteVolumeDistance(rightMSV);
}


static float dimension = 64; //yuck
float hausdorffDist(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV)
{
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leftMIV->hausdorffDistance(rightMIV, dimension);
	}

	return leftMIV->hausdorffDistance(rightMIV, dimension) + leftMSV->hausdorffDistance(rightMSV, dimension);
}

shapeMetricFn shapeDistance = averageRelVolume;

void setDistance(DistanceFunction df, float dim)
{
	switch(df)
	{
	case RelativeVolume:
		 shapeDistance = averageRelVolume;
		 break;
	case AbsVolume:
		 shapeDistance = averageAbsVolume;
		 break;
	case Hausdorff:
		dimension = dim;
		shapeDistance = hausdorffDist;
		break;
	case RelativeTriple:
		shapeDistance = tripleRelative;
		break;
	case AbsoluteTriple:
		shapeDistance = tripleAbsolute;
		break;
	}
}

//return a "distance" between obj and MIV/MSV; the lower the distance
//the  more likely a node should be searched
//min and max should bookend the ultimate leaf volume distances
float searchVolumeDist(const MappableOctTree* obj, const MappableOctTree* MIV,
		const MappableOctTree* MSV, float& min, float& max)
{
	min = 1 - obj->intersectVolume(MSV)/obj->unionVolume(MIV); //percent of shape already covered
	max = 1 - obj->intersectVolume(MIV)/obj->unionVolume(MSV); //difference if MIV/MSV after merge

	//selectivity?
	return min + max;
}

//a nearest neighbors distance, should match up with bounds of searchVolumeDist
float volumeDist(const MappableOctTree* x, const MappableOctTree* y)
{
	return 1 - x->intersectVolume(y)/x->unionVolume(y);
}

