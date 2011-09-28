/*
 * PtrPtrOctTree.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#include "PtrOctTree.h"

#include <cassert>

//recursive creation of oct tree
void PtrOctTree::create(OctNode *node, const Cube& cube,
		const vector<MolSphere>& mol)
{
	//does the mol overlap with this cube?
	vector<MolSphere> pruned;
	pruned.reserve(mol.size());
	for (unsigned i = 0, n = mol.size(); i < n; i++)
	{
		if (cube.intersectsSphere(mol[i]))
		{
			pruned.push_back(mol[i]);
		}
	}

	if (pruned.size() == 0)
	{
		//no overlap, all done
		node->type = Empty;
	}
	else if (cube.getDimension() <= resolution) //consider it full
	{
		node->type = Full;
	}
	else //subdivide into children
	{
		node->type = Children;
		unsigned fullcnt = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			node->children[i] = new OctNode();
			Cube newc = cube.getOctant(i);
			create(node->children[i], newc, pruned);
			if (node->children[i]->type == Full)
				fullcnt++;
		}

		//are all the children full? then truncate and mark node as full
		if (fullcnt == 8)
		{
			node->deleteChildren();
			node->type = Full;
		}
	}
}

PtrOctTree::OctNode::OctNode(const OctNode& rhs) :
	type(Empty)
{
	*this = rhs;
}

PtrOctTree::OctNode& PtrOctTree::OctNode::operator=(const PtrOctTree::OctNode& rhs)
{
	deleteChildren();
	type = rhs.type;

	memset(children, 0, sizeof(children));

	if (type == Children)
	{
		for (unsigned i = 0; i < 8; i++)
		{
			children[i] = new OctNode(*rhs.children[i]);
		}
	}

	return *this;
}

//recursively change this node to be its inverse
void PtrOctTree::OctNode::invert()
{
	if (type == Empty)
	{
		type = Full;
	}
	else if (type == Full)
	{
		type = Empty;
	}
	else //children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			children[i]->invert();
		}
	}
}

//recursively change this node to be intersected with rhs
bool PtrOctTree::OctNode::intersect(const OctNode *rhs)
{
	if (rhs->type == Empty)
	{
		deleteChildren();
		return true;
	}
	else if (rhs->type == Full || type == Empty)
	{
		//identity
		return false;
	}
	else if (type == Full)
	{
		//just copy
		*this = *rhs;
		return true;
	}
	else //both children
	{
		unsigned numEmpty = 0;
		bool changed = false;
		for (unsigned i = 0; i < 8; i++)
		{
			changed |= children[i]->intersect(rhs->children[i]);
			if (children[i]->type == Empty)
				numEmpty++;
		}

		//if all children empty, can truncate without loss of precision
		if (numEmpty == 8)
		{
			deleteChildren();
		}
		return changed;
	}
}

bool PtrOctTree::OctNode::unionWith(const OctNode *rhs)
{
	if (type == Full || rhs->type == Empty) //just identiy
	{
		return false;
	}
	else if (rhs->type == Full)
	{
		deleteChildren();
		type = Full;
		return true;
	}
	else if (type == Empty)
	{
		//same as rhs
		*this = *rhs;
		return true;
	}
	else //both have children
	{
		unsigned numFull = 0;
		bool changed = false;
		for (unsigned i = 0; i < 8; i++)
		{
			changed |= children[i]->unionWith(rhs->children[i]);
			if (children[i]->type == Full)
				numFull++;
		}

		//if all full, can truncate without loss of precision
		if (numFull == 8)
		{
			deleteChildren();
			type = Full;
		}
		return changed;
	}
}

//this computes the maximum included volume at the specified resolution
void PtrOctTree::OctNode::truncate(double resolution, float dim)
{
	if (type == Children)
	{
		if (dim <= resolution)
		{
			//make blank
			deleteChildren();
		}
		else //unchanged, process children
		{
			for (unsigned i = 0; i < 8; i++)
			{
				children[i]->truncate(resolution, dim / 2);
			}
		}
	}
}

//computes the minimum surrounding volume at the specified resolution
void PtrOctTree::OctNode::grow(double resolution, float dim)
{
	if (type == Children)
	{
		if (dim <= resolution)
		{
			//make full
			deleteChildren();
			type = Full;
		}
		else //unchanged, process children
		{
			for (unsigned i = 0; i < 8; i++)
			{
				children[i]->truncate(resolution, dim / 2);
			}
		}
	}
}

//recursive helper for volume
float PtrOctTree::OctNode::volume(float dim) const
{
	if (type == Empty)
		return 0;
	else if (type == Full)
		return dim * dim * dim;

	//has children
	float ret = 0;
	for (unsigned i = 0; i < 8; i++)
	{
		ret += children[i]->volume(dim / 2);
	}
	return ret;
}

//recursive helper for numleaves
unsigned PtrOctTree::OctNode::leaves() const
{
	if (type != Children)
		return 1;

	unsigned ret = 0;
	for (unsigned i = 0; i < 8; i++)
	{
		ret += children[i]->leaves();
	}
	return ret;
}

//return volume of intersection
float PtrOctTree::OctNode::intersectVolume(const OctNode *rhs, float dim) const
{
	if (rhs->type == Empty)
	{
		return 0;
	}
	else if (rhs->type == Full || type == Empty)
	{
		return volume(dim);
	}
	else if (type == Full)
	{
		return rhs->volume(dim);
	}
	else //both children
	{
		float vol = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			vol += children[i]->intersectVolume(rhs->children[i], dim / 2);
		}
		return vol;
	}
}

//return volume of union
float PtrOctTree::OctNode::unionVolume(const OctNode *rhs, float dim) const
{
	if (type == Full || rhs->type == Empty) //just identiy
	{
		return volume(dim);
	}
	else if (rhs->type == Full || type == Empty)
	{
		return rhs->volume(dim);
	}
	else //both have children
	{
		float vol = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			vol += children[i]->unionVolume(rhs->children[i], dim / 2);
		}
		return vol;
	}
}

void PtrOctTree::OctNode::write(ostream& out) const
{

	out.write((char*) &type, sizeof(type));

	if (type == Children)
	{
		//all children must exist
		for (unsigned i = 0; i < 8; i++)
		{
			assert(children[i] != NULL);
			children[i]->write(out);
		}
	}
}

void PtrOctTree::OctNode::read(istream& in)
{
	deleteChildren();
	in.read((char*) &type, sizeof(type));

	if (type == Children)
	{
		for (unsigned i = 0; i < 8; i++)
		{
			children[i] = new OctNode();
			children[i]->read(in);
		}
	}
}

//dump oct tree to a file in binary form
void PtrOctTree::write(ostream& out) const
{
	out.write((char*) &dimension, sizeof(dimension));
	out.write((char*) &resolution, sizeof(resolution));

	root->write(out);
}

//overwrite the current octree
void PtrOctTree::read(istream& in)
{
	in.read((char*) &dimension, sizeof(dimension));
	in.read((char*) &resolution, sizeof(resolution));

	root->read(in);
}

