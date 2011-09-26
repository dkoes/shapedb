/*
 * OctTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "OctTree.h"
#include <cassert>

//recursive creation of oct tree
void OctTree::create(OctNode *node, const Cube& cube,
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

OctTree::OctNode::OctNode(const OctNode& rhs) :
	type(Empty)
{
	*this = rhs;
}

OctTree::OctNode& OctTree::OctNode::operator=(const OctTree::OctNode& rhs)
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
void OctTree::OctNode::invert()
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
bool OctTree::OctNode::intersect(const OctNode *rhs)
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

bool OctTree::OctNode::unionWith(const OctNode *rhs)
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
void OctTree::OctNode::truncate(double resolution, float dim)
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
void OctTree::OctNode::grow(double resolution, float dim)
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
float OctTree::OctNode::volume(float dim) const
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
unsigned OctTree::OctNode::leaves() const
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
float OctTree::OctNode::intersectVolume(const OctNode *rhs, float dim) const
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
float OctTree::OctNode::unionVolume(const OctNode *rhs, float dim) const
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

void OctTree::OctNode::write(ostream& out) const
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

void OctTree::OctNode::read(istream& in)
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
void OctTree::write(ostream& out) const
{
	out.write((char*) &dimension, sizeof(dimension));
	out.write((char*) &resolution, sizeof(resolution));

	root->write(out);
}

//overwrite the current octree
void OctTree::read(istream& in)
{
	in.read((char*) &dimension, sizeof(dimension));
	in.read((char*) &resolution, sizeof(resolution));

	root->read(in);
}


//helper functions for helping me deal with the lack of internal nodes

//return the position after a tree of the given level
unsigned LinearOctTree::absorbTreeAtLevel(const vector<OctVal>& T, unsigned pos, unsigned level)
{
	if(T[pos].level == level)
		return pos+1;
	//otherwise need to read 8 trees of the next higher level
	for(unsigned i = 0; i < 8; i++)
	{
		pos = absorbTreeAtLevel(T, pos, level+1);
	}
	return pos;
}

//like absorb, but append trees that are read
unsigned LinearOctTree::appendTreeAtLevel(const vector<OctVal>& T, unsigned pos, unsigned level, vector<OctVal>& appendto)
{
	if(T[pos].level == level)
	{
		appendto.push_back(T[pos]);
		return pos+1;
	}
	//otherwise need to read 8 trees of the next higher level
	for(unsigned i = 0; i < 8; i++)
	{
		pos = appendTreeAtLevel(T, pos, level+1, appendto);
	}
	return pos;
}

//like above, but just compute volume
unsigned LinearOctTree::volumeOfTreeAtLevel(const vector<OctVal>& T, unsigned pos, unsigned level, float& vol) const
{
	if(T[pos].level == level)
	{
		if(T[pos].flag == Full)
			vol += levelVolumes[level];
		return pos+1;
	}
	//otherwise need to read 8 trees of the next higher level
	for(unsigned i = 0; i < 8; i++)
	{
		pos = volumeOfTreeAtLevel(T, pos, level+1, vol);
	}
	return pos;
}

//precompute volume of a cube at a given level
void LinearOctTree::setLevelVolumes()
{
	levelVolumes.clear();
	float dim = dimension;
	float res = resolution/2; //just to be safe
	if(res <= 0) return;
	while(dim >= res)
	{
		float vol = dim*dim*dim;
		levelVolumes.push_back(vol);
		dim /= 2;
	}
}

//find the position that points to the subtree specified by coord, or the highest leaf
unsigned LinearOctTree::traverseOctantCoord(const vector<OctVal>& T, const vector<unsigned>& coord)
{
	unsigned pos = 0;
	for(unsigned i = 0, n = coord.size(); i < n; i++)
	{
		//i is the level, coord[i] is what octant in the next level to look at
		if(T[pos].level == i)
		{
			return pos;
		}
		else
		{
			for(unsigned j = 0, m = coord[i]; j < m; j++)
			{
				pos = absorbTreeAtLevel(T, pos, i+1);
			}
		}
	}
	return pos;
}

//there are 256 possible bit patterns within a single octant, this will return
//the appropriate pattern for the specified coordinates
unsigned LinearOctTree::getOctantPattern(const vector<unsigned>& coord) const
{
	//get a cursor to the correct position
	if(tree.size() == 0)
		return 0;
	unsigned level = coord.size();
	unsigned pos = traverseOctantCoord(tree, coord);

	if(tree[pos].level <= level)
	{
		//request octant is subsummed
		if(tree[pos].flag == Empty)
			return 0;
		else
			return 255;
	}
	else
	{
		unsigned nextlevel = level+1;
		unsigned ret = 0;
		for(unsigned i = 0; i < 8; i++)
		{
			if(tree[pos].level == nextlevel) //exactly matches an octant
			{
				if(tree[pos].flag == Full)
					ret |= (1<<i);
			}
			else if(tree[pos].level > nextlevel) //something there
			{
				ret |= (1<<i);
			}
			pos = absorbTreeAtLevel(tree, pos, nextlevel);
		}
		return ret;
	}
}


//create a linear oct tree recursively
void LinearOctTree::create(const Cube& cube, unsigned level,
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
		tree.push_back(OctVal(Empty, level));
	}
	else if (cube.getDimension() <= resolution) //consider it full
	{
		tree.push_back(OctVal(Full, level));
	}
	else //subdivide into children
	{
		unsigned curend = tree.size();
		unsigned fullcnt = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			Cube newc = cube.getOctant(i);
			create(newc, level+1, pruned);
			if(tree.back().level == level+1 && tree.back().flag == Full)
				fullcnt++;
		}

		//are all the children full? then truncate and mark node as full
		if (fullcnt == 8)
		{
			assert(tree.size() - curend == 8);
			tree.resize(tree.size()-8);
			tree.push_back(OctVal(Full,level));
		}
	}
}

//the volume, can just traverse the whole tree
float LinearOctTree::volume() const
{
	float vol = 0;
	for(unsigned i = 0, n = tree.size(); i < n; i++)
	{
		if(tree[i].flag == Full)
		{
			vol += levelVolumes[tree[i].level];
		}
	}
	return vol;
}

//invert filled and unfilled
void LinearOctTree::invert()
{
	//can simply reverse flags since interior nodes
	//are unchanged in inverted tree
	for(unsigned i = 0, n = tree.size(); i < n; i++)
	{
		if(tree[i].flag == Empty)
			tree[i].flag = Full;
		else
			tree[i].flag = Empty;
	}
}



//mogrifying union or intersection - code is almost the same
bool LinearOctTree::operation(const LinearOctTree& rhs, bool doUnion)
{
	vector<OctVal> result;
	if(doUnion)
		result.reserve(max(rhs.tree.size(), tree.size()));
	else
		result.reserve(min(rhs.tree.size(), tree.size()));

	Type op = (doUnion ? Full : Empty);
	Type notop = (doUnion ? Empty : Full);

	unsigned lpos = 0, rpos = 0;
	bool changed = false;
	//simultaneously iterate over both trees
	assert(dimension == rhs.dimension && resolution == rhs.resolution);
	while(lpos < tree.size() && rpos < rhs.tree.size())
	{
		unsigned llevel = tree[lpos].level;
		unsigned rlevel = rhs.tree[rpos].level;
		Type ltype = tree[lpos].flag;
		Type rtype = rhs.tree[rpos].flag;

		if(llevel == rlevel) //same level
		{
			if(ltype == op || rtype == op)
			{
				result.push_back(OctVal(op, llevel));
				if(ltype != op)
					changed = true;
			}
			else //both empty
			{
				result.push_back(OctVal(notop, llevel));
			}
			lpos++;
			rpos++;
		}
		else if(llevel < rlevel)
		{
			if(ltype == op)
			{
				result.push_back(OctVal(op, llevel));
				rpos = absorbTreeAtLevel(rhs.tree, rpos, llevel);
			}
			else
			{
				//result equals rhs
				rpos = appendTreeAtLevel(rhs.tree, rpos, llevel, result);
				changed = true;
			}
			lpos++;
		}
		else if(llevel > rlevel)
		{
			if(rtype == op)
			{
				result.push_back(OctVal(op, rlevel));
				lpos = absorbTreeAtLevel(tree, lpos, rlevel);
				changed = true;
			}
			else
			{
				//result equals lhs
				lpos = appendTreeAtLevel(tree, lpos, rlevel, result);
			}
			rpos++;
		}
	}
	assert(rpos == rhs.tree.size());
	assert(lpos == tree.size());

	swap(tree, result);
	return changed;
}


//return the volume of the result of a union/intersect
float LinearOctTree::volOperation(const LinearOctTree& rhs, bool doUnion) const
{
	float vol = 0;
	Type op = (doUnion ? Full : Empty);
	Type notop = (doUnion ? Empty : Full);

	unsigned lpos = 0, rpos = 0;
	//simultaneously iterate over both trees
	assert(dimension == rhs.dimension && resolution == rhs.resolution);
	unsigned lsize = tree.size();
	unsigned rsize = rhs.tree.size();
	while(lpos < lsize && rpos < rsize)
	{
		unsigned llevel = tree[lpos].level;
		unsigned rlevel = rhs.tree[rpos].level;
		Type ltype = tree[lpos].flag;
		Type rtype = rhs.tree[rpos].flag;

		if(llevel == rlevel) //same level
		{
			if(ltype == op || rtype == op)
			{
				if(op == Full)
					vol += levelVolumes[llevel];
			}
			else //both empty or both full
			{
				if(notop == Full)
					vol += levelVolumes[llevel];
			}
			lpos++;
			rpos++;
		}
		else if(llevel < rlevel)
		{
			if(ltype == op)
			{
				if(op == Full)
					vol += levelVolumes[llevel];
				rpos = absorbTreeAtLevel(rhs.tree, rpos, llevel);
			}
			else
			{
				//result equals rhs
				rpos = volumeOfTreeAtLevel(rhs.tree, rpos, llevel, vol);
			}
			lpos++;
		}
		else if(llevel > rlevel)
		{
			if(rtype == op)
			{
				if(op == Full)
					vol += levelVolumes[rlevel];
				lpos = absorbTreeAtLevel(tree, lpos, rlevel);
			}
			else
			{
				//result equals lhs
				lpos = volumeOfTreeAtLevel(tree, lpos, rlevel, vol);
			}
			rpos++;
		}
	}
	assert(rpos == rhs.tree.size());
	assert(lpos == tree.size());

	return vol;
}


//mogrifying intersection
bool LinearOctTree::intersect(const LinearOctTree& rhs)
{
	if(&rhs == this)
		return false;
	return operation(rhs, false);
}

bool LinearOctTree::unionWith(const LinearOctTree& rhs)
{
	if(&rhs == this)
		return false;
	return operation(rhs, true);
}


//volume calculations that don't require creating a tmp tree
float LinearOctTree::intersectVolume(const LinearOctTree& rhs) const
{
	return volOperation(rhs, false);
}

float LinearOctTree::unionVolume(const LinearOctTree& rhs) const
{
	return volOperation(rhs, true);
}


void LinearOctTree::write(ostream& out) const
{
	out.write((char*)&dimension, sizeof(dimension));
	out.write((char*)&resolution, sizeof(resolution));
	unsigned n = tree.size();
	out.write((char*)&n, sizeof(n));
	for(unsigned i = 0; i < n; i++)
	{
		out.write((char*)&tree[i], sizeof(tree[i]));
	}
}

void LinearOctTree::read(istream& in)
{
	in.read((char*)&dimension, sizeof(dimension));
	in.read((char*)&resolution, sizeof(resolution));

	setLevelVolumes();

	unsigned n = 0;
	in.read((char*)&n, sizeof(n));

	tree.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		in.read((char*)&tree[i], sizeof(tree[i]));
	}
}
