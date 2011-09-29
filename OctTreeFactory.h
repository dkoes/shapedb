/*
 * OctTreeFactory.h
 *
 * Factory for creating different types of oct trees.  A user class
 * need not worry what kinds is generate, although the assumption is
 * that different kinds will not be mixed.
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#ifndef OCTTREEFACTORY_H_
#define OCTTREEFACTORY_H_

#include "OctTree.h"


class OctTreeFactory
{
public:
	enum OctTreeType {Linear, Pointer,Array};
private:
	OctTreeType type;
public:
	OctTreeFactory(): type(Linear) {}
	OctTreeFactory(OctTreeType t): type(t) {}


	OctTree* newOctTree() const { return newOctTree(0,0); }
	OctTree* newOctTree(float dim, float res) const {
		vector<MolSphere> dummy;
		return newOctTree(dim, res, dummy);
	}

	OctTree* newOctTree(float dim, float res,
			const vector<MolSphere>& mol) const;

	void write(ostream& out) const
	{
		out.write((char*)&type, sizeof(type));
	}

	void read(istream& in)
	{
		in.read((char*)&type, sizeof(type));
	}
};

#endif /* OCTTREEFACTORY_H_ */
