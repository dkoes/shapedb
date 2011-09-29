/*
 * OctTreeFactory.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#include "OctTreeFactory.h"
#include "LinearOctTree.h"
#include "PtrOctTree.h"
#include "ArrayOctTree.h"

OctTree* OctTreeFactory::newOctTree(float dim, float res,
		const vector<MolSphere>& mol) const
{
	switch (type)
	{
	case Linear:
		return new LinearOctTree(dim, res, mol);
	case Pointer:
		return new PtrOctTree(dim, res, mol);
	case Array:
		return new ArrayOctTree(dim, res, mol);
	}
	return NULL;
}
