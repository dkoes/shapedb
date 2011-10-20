/*
 * MappableOctTree.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "MappableOctTree.h"
#include <cstring>
#include <cassert>

MappableOctTree* MappableOctTree::clone() const
{
	unsigned sz = bytes();
	void* mem = malloc(sz);
	memcpy(mem, this, sz);
	return (MappableOctTree*)mem;
}



MChildNode MappableOctTree::createFrom_r(unsigned N, MChildNode* nodes, const MappableOctTree** trees, vector<MOctNode>& newtree, bool isUnion)
{
	float bitvol = 0;
	unsigned numFilled = 0;
	unsigned andpat = 0xff, orpat = 0;
	for(unsigned i = 0; i < N; i++)
	{
		if(nodes[i].isLeaf)
		{
			numFilled++;
			andpat &= nodes[i].pattern;
			orpat |= nodes[i].pattern;

			if(nodes[i].pattern != 0 && bitvol == 0)
			{
				//compute volume of a single bit pattern
				bitvol = nodes[i].volume()/nodes[i].index;
			}
		}
	}

	if(isUnion && (numFilled == N || orpat == 0xff)) //just return a patterned node
	{
		unsigned nb = __builtin_popcount(orpat);
		MChildNode ret(true, orpat,nb, bitvol*nb);
		return ret;
	}
	else if(!isUnion && (numFilled == N || andpat == 0))
	{
		unsigned nb = __builtin_popcount(andpat);
		MChildNode ret(true, andpat,nb, bitvol*nb);
		return ret;
	}
	else //need to look at non-leaf children
	{
		unsigned pos = newtree.size();
		MChildNode ret(false, 0, pos, 0);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		MChildNode nextnodes[N+1];
		const MappableOctTree *nexttrees[N+1];

		//setup nonleaf trees
		unsigned n = 0;
		for(unsigned i = 0; i < N; i++)
		{
			if(!nodes[i].isLeaf)
			{
				nexttrees[n] = trees[i];
				n++;
			}
		}
		if(numFilled > 0)
			nexttrees[n++] = NULL; //dummy for merge of leaves
		for (unsigned i = 0; i < 8; i++)
		{
			//create array of sub-octants from non-leafs
			unsigned cnt = 0;
			for(unsigned j = 0; j < N; j++)
			{
				if(!nodes[j].isLeaf)
				{
					nextnodes[cnt] = trees[j]->tree[nodes[j].index].children[i];
					cnt++;
				}
			}
			//summerize filled nodes
			if(numFilled > 0)
			{
				unsigned pat = isUnion ? orpat : andpat;
				if(pat & (1<<i))
					pat = 0xff;
				else
					pat = 0;
				nextnodes[cnt] = MChildNode(true, pat, pat ? 8 : 0, pat ? bitvol : 0);
				cnt++;
			}
			assert(cnt == n);

			MChildNode nchild = createFrom_r(n, nextnodes, nexttrees, newtree, isUnion);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && (nchild.pattern == 0 || nchild.pattern == 0xff))
			{
				numFullEmpty++;
				if(nchild.pattern == 0xff)
					newPattern |= 1<<i;
			}

			ret.vol += nchild.volume();
		}

		if(numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.pattern = newPattern;
			ret.index = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

//invert tree, this can be done inplace
void MappableOctTree::invert(float dim)
{
	float expectedV = dim*dim*dim-root.volume();
	root.invert(tree, dim*dim*dim);
	assert(root.volume() == expectedV);
}

//union
MappableOctTree* MappableOctTree::createFromUnion(unsigned N, const MappableOctTree** in)
{
	MChildNode roots[N];
	for(unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in,  newtree, true);

	unsigned sz = sizeof(MappableOctTree)+newtree.size()*sizeof(MOctNode);
	void *mem = malloc(sz);
	return new (mem) MappableOctTree(newroot, newtree);
}

MappableOctTree* MappableOctTree::createFromIntersection(unsigned N, const MappableOctTree** in)
{
	MChildNode roots[N];
	for(unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in,  newtree, false);

	unsigned sz = sizeof(MappableOctTree)+newtree.size()*sizeof(MOctNode);
	void *mem = malloc(sz);
	return new (mem) MappableOctTree(newroot, newtree);
}


//volume calculations that don't require creating a tmp tree
float MappableOctTree::intersectVolume(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);
	return ival;
}

float MappableOctTree::unionVolume(const MappableOctTree *rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);
	return uval;
}


bool MappableOctTree::containedIn(const MappableOctTree *rhs) const
{
	return root.containedIn(tree, rhs->root, rhs->tree);
}

//return total volume contained in octtree
float MappableOctTree::volume() const
{
	return root.volume();
}

//return number of leaves
unsigned MappableOctTree::leaves() const
{
	unsigned cnt = 0;
	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < 8; j++)
		{
			if (tree[i].children[j].isLeaf)
				cnt++;
		}
	}
	return cnt;
}


void MappableOctTree::write(ostream& out) const
{
	unsigned sz = sizeof(MappableOctTree)+N*sizeof(MOctNode);
	out.write((char*)this, sz);
}

float MappableOctTree::relativeVolumeDistance(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);

	return 1 - ival/uval;
}

float MappableOctTree::absoluteVolumeDistance(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree, rhs->root,rhs->tree, ival, uval);

	return uval - ival;
}



void MChildNode::invert(MOctNode* tree, float maxvol)
{
	if(isLeaf)
	{
		pattern = ~pattern;
		index = 8 - index;
		vol = maxvol - vol;
	}
	else
	{
		vol = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].invert(tree, maxvol/8);
			vol += tree[index].children[i].vol;
		}
	}
}


//compute both the intersection and union volume at once
void MChildNode::intersectUnionVolume(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree, float& intersectval, float& unionval) const

{
	if(rhs.isLeaf && isLeaf)
	{
		if(pattern == 0 || rhs.pattern == 0)
			intersectval += 0; //no intersect val
		else
		{
			float bitv = vol/index; //volume of each bit
			unsigned intpat = pattern & rhs.pattern;

			intersectval += __builtin_popcount(intpat)*bitv;
		}
		if(pattern == 0)
			unionval += rhs.volume();
		else if(rhs.pattern == 0)
			unionval += volume();
		else
		{
			float bitv = vol/index; //volume of each bit
			unsigned orpat = pattern | rhs.pattern;

			unionval += __builtin_popcount(orpat)*bitv;
		}
	}
	else if(isLeaf && pattern == 0xff)
	{
		unionval += volume();
		intersectval += rhs.volume();
	}
	else if (rhs.isLeaf && rhs.pattern == 0xff)
	{
		unionval += rhs.volume();
		intersectval += volume();
	}
	else if (rhs.isLeaf && rhs.pattern == 0)
	{
		//vol of this tree
		unionval += volume();
	}
	else if (isLeaf && pattern == 0)
	{
		// from rhs
		unionval += rhs.volume();
	}
	else if(isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(index > 0);
		//rhs has children, have to compare bit by bit
		float bvol = vol/index;
		for(unsigned i = 0; i < 8; i++)
		{
			MChildNode rchild = rtree[rhs.index].children[i];
			if(pattern & (1<<i))
			{
				unionval += bvol;
				intersectval += rchild.volume();
			}
			else
			{
				unionval += rchild.volume();
			}
		}
	}
	else if(rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.index > 0);
		//rhs has children, have to compare bit by bit
		float bvol = rhs.vol/rhs.index;
		for(unsigned i = 0; i < 8; i++)
		{
			MChildNode child = tree[index].children[i];
			if(rhs.pattern & (1<<i))
			{
				unionval += bvol;
				intersectval += child.volume();
			}
			else
			{
				unionval += child.volume();
			}
		}
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].intersectUnionVolume(tree,
					rtree[rhs.index].children[i], rtree, intersectval,unionval);
		}
	}
}

//return true if this is contained in rhs - short circuit eval
bool MChildNode::containedIn(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const
{
	if(rhs.isLeaf && isLeaf)
	{
		return (pattern & rhs.pattern) == pattern;
	}
	else if(isLeaf && pattern == 0)
	{
		return true;
	}
	else if (rhs.isLeaf && rhs.pattern == 0)
	{
		return false; //somthing not in nothing (nothing handled above)
	}
	else if (rhs.isLeaf && rhs.pattern == 0xff)
	{
		return true;
	}
	else if (isLeaf && pattern == 0xff)
	{
		return false;
	}
	else if(isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(index > 0);
		//rhs has children, have to compare bit by bit
		for(unsigned i = 0; i < 8; i++)
		{
			if(pattern & (1<<i))
			{
				MChildNode rchild = rtree[rhs.index].children[i];
				//rchild must be a leaf and full
				if(!rchild.isLeaf || rchild.pattern != 0xff)
					return false;
			}
		}
	}
	else if(rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.index > 0);
		//rhs has children, have to compare bit by bit
		for(unsigned i = 0; i < 8; i++)
		{
			if(!(rhs.pattern & (1<<i)))
			{
				MChildNode child = tree[index].children[i];
				//if rhs's bit is empty, then our child must be empty
				if(!child.isLeaf || child.pattern != 0)
					return false;
			}
		}
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if(!tree[index].children[i].containedIn(tree,
					rtree[rhs.index].children[i], rtree))
				return false;
		}
		return true;;
	}
	return false;
}

