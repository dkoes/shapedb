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
float (*shapeDistanceFn)(const MappableOctTree* leftMIV, const MappableOctTree* leftMSV,
		const MappableOctTree* rightMIV, const MappableOctTree* rightMSV);

//a global function for performing shape comparisons between MIV/MSV pairs
extern shapeDistanceFn shapeDistance;

#endif /* SHAPEDISTANCE_H_ */
