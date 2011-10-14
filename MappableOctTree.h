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

#include "ArrayOctTree.h"

struct MOctNode;
struct MChildNode
{
	bool isLeaf :1;
	bool isFull :1; //only relevant if leaf
	unsigned index :30;
	float vol;

	MChildNode() :
		isLeaf(true), isFull(false), index(0), vol(-1)
	{
	}

	MChildNode(bool leaf, bool full) :
		isLeaf(leaf), isFull(full), index(0), vol(-1)
	{
	}
	MChildNode(bool leaf, bool full, unsigned i) :
		isLeaf(leaf), isFull(full), index(i), vol(-1)
	{
	}



	void intersectUnionVolume(const MOctNode* tree,
			const MOctNode* rtree, const MChildNode& rhs, float dim,
			float& intersectval, float& unionval) const;

	float
	intersectVolume(const MOctNode* tree, const MOctNode* rtree,
			const MChildNode& rhs, float dim) const;
	float
	unionVolume(const MOctNode* tree, const MOctNode* rtree,
			const MChildNode& rhs, float dim) const;

	bool containedIn(const MOctNode* tree, const MOctNode* rtree,
			const MChildNode& rhs) const;

	float
			subtractVolume(const MOctNode* tree,
					const MOctNode* bigtree, const MChildNode& bigrhs,
					float dim) const;

	float volume() const { return vol; }

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

	MappableOctTree(const MappableOctTree& rhs);

	static MappableOctTree* create_r(float dim, float res, const Cube& cube,
			const vector<MolSphere>& mol);

public:

	template <class Object>
	static MappableOctTree* create(float dim, float res, const Object& obj);
	static MappableOctTree* create(const ArrayOctTree& src);

	MappableOctTree* clone() const;

	//invert filled and unfilled - change in place
	void invert();

	//return intersect of the n trees found in arr
	MappableOctTree* createFromIntersection(const vector<const MappableOctTree*>& in);
	//union
	MappableOctTree* createFromUnion(const vector<const MappableOctTree*>& in);

	//volume calculations that don't require creating a tmp tree
	float intersectVolume(const MappableOctTree * rhs) const;
	float unionVolume(const MappableOctTree *rhs) const;
	float subtractVolume(const MappableOctTree *rhs) const;

	bool containedIn(const MappableOctTree *rhs) const;

	//return total volume contained in octtree
	float volume() const;

	//return number of leaves
	unsigned leaves() const;

	void write(ostream& out) const;

	float hausdorffDistance(const OctTree* B) const;
	float relativeVolumeDistance(const OctTree * B) const;
	float absoluteVolumeDistance(const OctTree * B) const;

};

#endif /* MAPPABLEOCTTREE_H_ */
