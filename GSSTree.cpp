/*
 * GSSTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "GSSTree.h"
#include "OctTree.h"
#include <math.h>
#include <iostream>


const unsigned GSSTree::MaxSplit = 8; //number of children in each node


void GSSTree::setBoundingBox(array<float,6>& box)
{
	float dims[3] = {0,0,0};
	for(unsigned i = 0; i < 3; i++)
	{
		min[i] = box[2*i];
		dims[i] = box[2*i+1]-min[i];
	}

	//to make things easier, cubify using the longest dimension
	dim = max(dims[2],max(dims[0],dims[1]));
	//round all dimensions up to a power of 2; use fp
	dim = pow(2, ceil(log(dim)/log(2)));
	if(dim < 1)
		dim = 1;
}

//add a single mol
void GSSTree::add(const vector<MolSphere>& m)
{
	//reposition mol
	vector<MolSphere> mol;
	mol.reserve(m.size());

	for(unsigned i = 0, n = m.size(); i < n; i++)
	{
		mol.push_back(MolSphere(m[i].x-min[0], m[i].y-min[1], m[i].z-min[2],m[i].r));
	}
	OctTree oct(dim, maxres, mol);

	//insert into the tree - recursively traverse full nodes along the most similar
	//path, the insert into a partially full node.  If no such node exists, split
	//the bottommost full node

	cerr << oct.volume() << " " << oct.leaves() << "\n";
}

//nearest neighbor search, return closest set of molspheres
void GSSTree::nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res)
{
	cerr << "Unimplemented: nn_search\n";
	exit(-1);
}


void GSSTree::write(const filesystem::path& out)
{
	cerr << "Unimplemented: write\n";
	exit(-1);}
