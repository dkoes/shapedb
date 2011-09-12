/*
 * GSSTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "GSSTree.h"
#include <math.h>


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
void GSSTree::add(const vector<MolSphere>& mol)
{

}

//nearest neighbor search, return closest set of molspheres
void GSSTree::nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res)
{
	abort(); //TODO
}


void GSSTree::write(const filesystem::path& out)
{
	abort(); //TODO
}
