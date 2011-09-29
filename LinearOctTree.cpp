/*
 * LinearOctTree.cpp
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#include "LinearOctTree.h"
#include <cassert>

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
unsigned LinearOctTree::getOctantPattern(const vector<unsigned>& coord, bool MSV) const
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
				if(MSV)
					ret |= (1<<i); //has some full children, but not all
			}
			pos = absorbTreeAtLevel(tree, pos, nextlevel);
		}
		return ret;
	}
}

//count the number of nodes at the specified level that have children
//(haven't been subsumed by larger nodes and aren't all full or all empty)
unsigned LinearOctTree::countInteriorNodesAtLevel(unsigned level) const
{
	unsigned ret = 0;
	unsigned pos = 0;
	while(pos < tree.size())
	{
		if(tree[pos].level == level)
			pos++;
		else if(tree[pos].level < level)
			pos++;
		else
		{
			//has children, absorb
			pos = absorbTreeAtLevel(tree, pos, level);
			ret++;
		}
	}
	return ret;
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
bool LinearOctTree::intersect(const OctTree* Rhs)
{
	const LinearOctTree& rhs = dynamic_cast<const LinearOctTree&>(*Rhs);
	if(&rhs == this)
		return false;
	return operation(rhs, false);
}

bool LinearOctTree::unionWith(const OctTree* Rhs)
{
	const LinearOctTree& rhs = dynamic_cast<const LinearOctTree&>(*Rhs);
	if(&rhs == this)
		return false;
	return operation(rhs, true);
}


//volume calculations that don't require creating a tmp tree
float LinearOctTree::intersectVolume(const OctTree* Rhs) const
{
	const LinearOctTree& rhs = dynamic_cast<const LinearOctTree&>(*Rhs);
	return volOperation(rhs, false);
}

float LinearOctTree::unionVolume(const OctTree* Rhs) const
{
	const LinearOctTree& rhs = dynamic_cast<const LinearOctTree&>(*Rhs);
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

//compute the directional Hausdorff distance,
//the maximum of the minimum distances from this to B
float LinearOctTree::hausdorffDistance(const OctTree* Bptr) const
{
	const LinearOctTree& B = dynamic_cast<const LinearOctTree&>(*Bptr);

	//for this oct tree representation, not sure we can get more efficient
	//than n^2 (easily)- read all the leaves of this tree and find the distance
	//to each leaf of B, finding the min

	//TODO: BUG: this is not actually accurate since we do not take
	//the context of the leaf cubes into account (only outside vertices should matter)


	float max = 0;
	for(LeafCubeIterator a(*this); a; ++a)
	{
		float min = HUGE_VAL;
		Cube aC = *a;
		for(LeafCubeIterator b(B); b; ++b)
		{
			Cube bC = *b;

			float d = aC.minDist(bC);
			if(d < min)
				min = d;
		}
		if(min > max)
			max = min;
	}
	return max;
}

//progress to the next cube
void LinearOctTree::LeafCubeIterator::step()
{
	do
	{
		pos++; //next leaf
		//remove current cube
		nestedCubes.pop_back();

		while(nestedOctants.back() == 7)
		{
			curlevel--; //zoom out
			nestedCubes.pop_back();
			nestedOctants.pop_back();
		}

		if(pos >= tree.size())
			break;

		nestedOctants.back()++;
		nestedCubes.push_back(nestedCubes.back().getOctant(nestedOctants.back()));


		assert(curlevel <= tree[pos].level);

		while(curlevel < tree[pos].level)
		{
			curlevel++;
			nestedCubes.push_back(nestedCubes.back().getOctant(0));
			nestedOctants.push_back(0);
		}

	} while(tree[pos].flag == Empty);

}
