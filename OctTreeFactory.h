/*
 * OctTreeFactory.h
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#ifndef OCTTREEFACTORY_H_
#define OCTTREEFACTORY_H_

#include "OctTree.h"
#include "LinearOctTree.h"

class OctTreeFactory
{
public:
	OctTreeFactory();
	virtual ~OctTreeFactory();

	OctTree* newOctTree() const { return newOctTree(0,0); }
	OctTree* newOctTree(float dim, float res) const {
		vector<MolSphere> dummy;
		return newOctTree(dim, res, dummy);
	}

	OctTree* newOctTree(float dim, float res, const vector<MolSphere>& mol) const
	{
		return new LinearOctTree(dim, res, mol);
	}

};

#endif /* OCTTREEFACTORY_H_ */
