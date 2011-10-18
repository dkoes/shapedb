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
	unsigned numEmpty = 0, numFull = 0;


	float fullVol = 0;
	for(unsigned i = 0; i < N; i++)
	{
		if(nodes[i].isLeaf)
		{
			if(nodes[i].isFull)
			{
				numFull++;
				fullVol = nodes[i].volume(); //this should be the same for all at this level
			}
			else
				numEmpty++;
		}
	}

	if((isUnion && numFull > 0) || numFull == N) //just return a full node
	{
		MChildNode ret(true, true, 0, fullVol);
		return ret;
	}
	else if((!isUnion && numEmpty > 0) || numEmpty == N) //return empty node
	{
		MChildNode ret(true, false, 0, 0);
		return ret;
	}
	else //need to look at non-leaf children
	{
		unsigned pos = newtree.size();
		MChildNode ret(false, false, pos, 0);
		newtree.push_back(MOctNode());

		unsigned numFull = 0;
		unsigned numEmpty = 0;
		MChildNode nextnodes[N];
		const MappableOctTree *nexttrees[N];

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
			assert(cnt == n);

			MChildNode nchild = createFrom_r(n, nextnodes, nexttrees, newtree, isUnion);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf)
			{
				if(nchild.isFull) numFull++;
				else numEmpty++;
			}

			ret.vol += nchild.volume();
		}

		if(numEmpty == n)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.isFull = false;
			ret.index = 0;
			ret.vol = 0;
		}
		else if(numFull == n)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.isFull = true;
			ret.index = 0;
			//ret.vol is correct
		}
		return ret;
	}
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
	return root.intersectVolume(tree, rhs->root, rhs->tree);
}

float MappableOctTree::unionVolume(const MappableOctTree *rhs) const
{
	return root.unionVolume(tree, rhs->root, rhs->tree);
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

//compute the interesection volume
float MChildNode::intersectVolume(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
			return volume();
		else
			return 0;
	}
	else if(isLeaf && !isFull)
	{
		return 0;
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		return 0;
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		//same as now
		return volume();
	}
	else if (isLeaf && isFull)
	{
		//same as rhs
		return rhs.volume();
	}
	else //both have children
	{
		float ret = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			ret += tree[index].children[i].intersectVolume(tree,
					rtree[rhs.index].children[i], rtree);
		}

		return ret;
	}
}

float MChildNode::unionVolume(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
			return volume();
		else
			return 0;
	}
	else if(isLeaf && isFull)
	{
		return volume();
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		return rhs.volume();
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		//vol of this tree
		return volume();
	}
	else if (isLeaf && !isFull)
	{
		// from rhs
		return rhs.volume();
	}
	else //both have children
	{
		float ret = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			ret += tree[index].children[i].unionVolume(tree,
					rtree[rhs.index].children[i], rtree);
		}

		return ret;
	}
}



//compute both the intersection and union volume at once
void MChildNode::intersectUnionVolume(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree, float& intersectval, float& unionval) const

{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
		{
			intersectval += volume();
			unionval += volume();
		}
	}
	else if(isLeaf && isFull)
	{
		unionval += volume();
		intersectval += rhs.volume();
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		unionval += rhs.volume();
		intersectval += volume();
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		//vol of this tree
		unionval += volume();
	}
	else if (isLeaf && !isFull)
	{
		// from rhs
		unionval += rhs.volume();
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
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		return true;
	}
	else if(isLeaf && !isFull)
	{
		return true;
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		return false; //somthing not in nothing (nothing handled above)
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		return true;
	}
	else if (isLeaf && isFull)
	{
		return false;
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

