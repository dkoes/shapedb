/*
 * MappableOctTree.h
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 *
 *  Based off of the array oct tree (several variants were tried) this
 *  is a mappable version (no virtual functions, no pointers)
 */

#ifndef MAPPABLEOCTTREE_H_
#define MAPPABLEOCTTREE_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "Cube.h"

using namespace std;

struct MOctNode;
struct MChildNode
{
	bool isLeaf :1;
	unsigned pattern: 8; //leaf pattern
	unsigned index :23; //ptr to next node, or bit cnt if leaf
	float vol; //this doubles the size of the tree, but makes similarity computation much faster, is the trade off worth it?

	MChildNode() :
			isLeaf(true), pattern(0), index(0), vol(-1)
	{
	}


	MChildNode(bool leaf, unsigned pat, unsigned i, float v) :
			isLeaf(leaf), pattern(pat), index(i), vol(v)
	{
	}
	void intersectUnionVolume(const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree,
			 float& intersectval,
			float& unionval) const;

	float
	intersectVolume(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;
	float
	unionVolume(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;

	void invert(MOctNode* tree, float maxvol);

	bool containedIn(const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const;

	float volume() const
	{
		return vol;
	}

};

struct MOctNode
{
	MChildNode children[8];

	MOctNode()
	{
	}

	MOctNode(const MChildNode& init)
	{
		std::fill_n(children, 8, init);
	}
};

class MappableOctTree
{
	MChildNode root;
	unsigned N; //number of octnodes
	MOctNode tree[]; //must memalloc beyond

	MappableOctTree()
	{
		//creation must be performed externally to properly allocate memory
	}

	//assume memory is allocate to tree
	MappableOctTree(const MChildNode& r, const vector<MOctNode>& nodes) :
			root(r), N(nodes.size())
	{
		for (unsigned i = 0, n = nodes.size(); i < n; i++)
		{
			tree[i] = nodes[i];
		}
	}

	MappableOctTree(const MappableOctTree& rhs);

	static MChildNode createFrom_r(unsigned N, MChildNode* nodes, const MappableOctTree** trees, vector<MOctNode>& newtree, bool isUnion);

public:

	MappableOctTree* clone() const;

	void invert(float dim);

	//return intersect of the n trees found in arr
	static MappableOctTree* createFromIntersection(unsigned N, const MappableOctTree** in);
	//union
	static MappableOctTree* createFromUnion(unsigned N, const MappableOctTree** in);

	//volume calculations that don't require creating a tmp tree
	float intersectVolume(const MappableOctTree * rhs) const;
	float unionVolume(const MappableOctTree *rhs) const;

	bool containedIn(const MappableOctTree *rhs) const;

	//return total volume contained in octtree
	float volume() const;

	//return number of leaves
	unsigned leaves() const;

	void write(ostream& out) const;

	float hausdorffDistance(const MappableOctTree* B) const;
	float relativeVolumeDistance(const MappableOctTree * B) const;
	float absoluteVolumeDistance(const MappableOctTree * B) const;

	unsigned bytes() const
	{
		return sizeof(MappableOctTree) + N*sizeof(MOctNode);
	}

private:
	template<class Object>
	static MChildNode create_r(float res, const Cube& cube, const Object& obj,
			vector<MOctNode>& tree)
	{
		//does the object overlap with this cube?
		bool intersects = obj.intersects(cube);

		MChildNode ret;
		if (!intersects)
		{
			//no overlap, all done
			ret.isLeaf = true;
			ret.pattern = 0;
			ret.index = 0;
			ret.vol = 0;
		}
		else if (cube.getDimension() <= res) //consider it full
		{
			ret.isLeaf = true;
			ret.pattern = 0xff;
			ret.index = 8;
			ret.vol = cube.volume();
		}
		else //subdivide into children
		{
			//assume this is an interior node
			ret.isLeaf = false;
			unsigned filledcnt = 0;
			unsigned pat = 0;
			unsigned bitcnt = 0;
			//allocate space
			unsigned pos = tree.size();
			ret.index = tree.size();
			tree.push_back(MOctNode());
			ret.vol = 0;
			for (unsigned i = 0; i < 8; i++)
			{
				Cube newc = cube.getOctant(i);
				MChildNode child = create_r(res, newc, obj, tree);
				tree[pos].children[i] = child;
				if (child.isLeaf)
				{
					if(child.pattern == 0)
						filledcnt++;
					else if(child.pattern == 0xff)
					{
						filledcnt++;
						pat |= (1<<i);
						bitcnt++;
					}
				}
				ret.vol += child.volume();
			}

			//are all the children full or empty? then can represent as a pattern
			if (filledcnt == 8)
			{
				tree.pop_back();
				ret.isLeaf = true;
				ret.pattern = pat;
				//volume is correct
				ret.index = bitcnt;
			}
		}
		return ret;
	}

public:

	//return a mappable oct tree, the pointer should be deallocated using free
	template<class Object>
	static MappableOctTree* create(float dim, float res, const Object& obj)
	{
		vector<MOctNode> nodes;
		MChildNode root = create_r(res, Cube(dim, -dim / 2, -dim / 2, -dim / 2),
				obj, nodes);
		unsigned N = nodes.size();
		void *mem = malloc(sizeof(MappableOctTree) + N * sizeof(MOctNode));

		MappableOctTree *ret = new (mem) MappableOctTree(root, nodes);

		return ret;
	}
};

#endif /* MAPPABLEOCTTREE_H_ */
