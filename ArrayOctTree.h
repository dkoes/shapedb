/*
 * ArrayOctTree.h
 *
 * This is a fully represented oct tree with (only!) internal nodes and "pointer"
 * but all the nodes are stored in a single vector and the pointers are indices
 *
 *  Created on: Sep 29, 2011
 *      Author: dkoes
 */

#ifndef ARRAYOCTTREE_H_
#define ARRAYOCTTREE_H_

#include "OctTree.h"

class ArrayOctTree: public OctTree
{
	float dimension; //sits within a cube of this dimension, power of 2
	float resolution; //increase resolution down to this increment (power of 2)

	enum OctType
	{
		Empty, Full, Interior
	};

	struct OctNode;
	struct ChildNode
	{
		bool isLeaf :1;
		bool isFull :1; //only relevant if leaf
		unsigned index :30;

		ChildNode() :
			isLeaf(true), isFull(false), index(0)
		{
		}

		ChildNode(bool leaf, bool full) :
			isLeaf(leaf), isFull(full), index(0)
		{
		}
		ChildNode(bool leaf, bool full, unsigned i) :
			isLeaf(leaf), isFull(full), index(i)
		{
		}

		ChildNode intersect(const vector<OctNode>& tree,
				const vector<OctNode>& rtree, const ChildNode& rhs,
				vector<OctNode>& newtree, bool& changed) const;
		ChildNode unionWith(const vector<OctNode>& tree,
				const vector<OctNode>& rtree, const ChildNode& rhs,
				vector<OctNode>& newtree, bool& changed) const;

		float
				intersectVolume(const vector<OctNode>& tree,
						const vector<OctNode>& rtree, const ChildNode& rhs,
						float dim) const;
		float
				unionVolume(const vector<OctNode>& tree,
						const vector<OctNode>& rtree, const ChildNode& rhs,
						float dim) const;

		bool containedIn(const vector<OctNode>& tree,
				const vector<OctNode>& rtree, const ChildNode& rhs) const;

		float volume(const vector<OctNode>& tree, float dim) const;
		unsigned getBitPattern(const vector<OctNode>& tree, bool MSV) const;

		ChildNode copyTo(const vector<OctNode>& from, vector<OctNode>& to) const;
	};

	struct OctNode
	{
		ChildNode children[8];

		OctNode()
		{
		}

		OctNode(const ChildNode& init)
		{
			std::fill_n(children, 8, init);
		}
	};

	vector<OctNode> tree;
	ChildNode root;
	ChildNode create(const Cube& cube, const vector<MolSphere>& mol);

public:
	ArrayOctTree()
	{
	}
	virtual ~ArrayOctTree()
	{
	}

	ArrayOctTree(const ArrayOctTree& rhs) :
		dimension(rhs.dimension), resolution(rhs.resolution)
	{

	}

	ArrayOctTree(float dim, float res, const vector<MolSphere>& mol) :
		dimension(dim), resolution(res)
	{
		Cube cube(dimension);
		root = create(cube, mol);
	}

	virtual OctTree* clone() const
	{
		return new ArrayOctTree(*this);
	}

	ArrayOctTree& operator=(ArrayOctTree rhs)
	{
		swap(*this, rhs); //rhs passed by value
		return *this;
	}

	friend void swap(ArrayOctTree& first, ArrayOctTree& second)
	{
		// enable ADL (not necessary in our case, but good practice)
		using std::swap;

		swap(first.resolution, second.resolution);
		swap(first.dimension, second.dimension);
		swap(first.tree, second.tree);
		swap(first.root, second.root);
	}

	//invert filled and unfilled
	virtual void invert();

	//mogrifying intersection
	virtual bool intersect(const OctTree* rhs);
	//mogrifying union
	virtual bool unionWith(const OctTree* rhs);

	//volume calculations that don't require creating a tmp tree
	virtual float intersectVolume(const OctTree * rhs) const;
	virtual float unionVolume(const OctTree *rhs) const;

	virtual bool containedIn(const OctTree *rhs) const;

	//return total volume contained in octtree
	virtual float volume() const;

	//return number of leaves
	virtual unsigned leaves() const;

	virtual void clear();
	virtual void fill();

	virtual void write(ostream& out) const;
	virtual void read(istream& in);

	virtual unsigned
			getOctantPattern(const vector<unsigned>& coord, bool MSV) const;

	virtual float hausdorffDistance(const OctTree* B) const;
};

#endif /* ARRAYOCTTREE_H_ */
