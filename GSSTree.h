/*
 * GSSTree.h
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 *
 *      http://dx.doi.org/10.1145/304181.304219
 *
 *      A GSS tree for molecular data.  When building takes as input
 *      a vector of spheres (x,y,z,r) that represent a molecule
 *      (eventaully will add meta data).
 */

#ifndef GSSTREE_H_
#define GSSTREE_H_

#include <boost/filesystem.hpp>
#include <boost/array.hpp>
#include <vector>

#include "MolSphere.h"

using namespace boost;
using namespace std;


//the gss tree
class GSSTree
{
	float maxres;
	float min[3];
	float dim;

	void setBoundingBox(array<float,6>& box);

public:
	GSSTree(filesystem::path& in); //initialize into memory from file; may mmap

	//set the global min and max extents of the GSS tree and
	//the highest resolution (default 1A)
	GSSTree(array<float,6>& boundingbox, double maxr=1.0): maxres(maxr)
	{
		setBoundingBox(boundingbox);
	}

	virtual ~GSSTree() {}

	//add a single mol
	void add(const vector<MolSphere>& mol);

	//nearest neighbor search, return closest set of molspheres
	void nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res);

	void write(const filesystem::path& out); //dump to file
};

#endif /* GSSTREE_H_ */
