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
extern shapeMetricFn shapeSimilarity;

#endif /* SHAPEDISTANCE_H_ */
