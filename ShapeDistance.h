/*
 * ShapeDistance.h
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 */

#ifndef SHAPEDISTANCE_H_
#define SHAPEDISTANCE_H_
#include "MappableOctTree.h"

typedef
float (*shapeMetricFn)(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV);

//a global function for performing shape comparisons between MIV/MSV pairs
extern shapeMetricFn shapeDistance;

enum DistanceFunction {RelativeVolume, AbsVolume, Hausdorff,RelativeTriple, AbsoluteTriple};

//set shapedistance to the requested function, hausdoff needs the dim
extern void setDistance(DistanceFunction df, float dim);

#endif /* SHAPEDISTANCE_H_ */
