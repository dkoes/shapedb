/*
 * PtrPtrOctTree.h
 *
 * This is a simple pointer-based oct tree where every octant is a standalone node.
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#ifndef PTROCTTREE_H_
#define PTROCTTREE_H_
#include "OctTree.h"

/* Pointer based oct tree. Tree is represented by node pointers
 */
class PtrOctTree
{
	float dimension; //sits within a cube of this dimension, power of 2
	float resolution; //increase resolution down to this increment (power of 2)

	enum OctType
	{
		Empty, Full, Children
	};
	struct OctNode
	{
		OctType type;
		OctNode *children[8]; //pointers belong to this node

		OctNode() :
			type(Empty)
		{
			memset(children, 0, sizeof(children));
		}

		OctNode(const OctNode& rhs);

		OctNode& operator=(const OctNode& rhs);

		void deleteChildren()
		{
			if (type == Children)
			{
				for (unsigned i = 0; i < 8; i++)
				{
					if (children[i])
					{
						delete children[i];
						children[i] = NULL;
					}
				}
			}
			type = Empty;
		}

		~OctNode()
		{
			deleteChildren();
		}

		//return true if changed
		void invert();
		bool intersect(const OctNode *rhs);
		bool unionWith(const OctNode *rhs);
		void truncate(double resolution, float dim);
		void grow(double resolution, float dim);

		float intersectVolume(const OctNode *rhs, float dim) const;
		float unionVolume(const OctNode *rhs, float dim) const;
		float volume(float dim) const; //recursive volume calculation
		unsigned leaves() const; //recursive leafs cnt

		void write(ostream& out) const;
		void read(istream& out);

	};

	OctNode *root;

	//recursively generate an octtree
	void create(OctNode *node, const Cube& cube, const vector<MolSphere>& mol);

public:

	PtrOctTree() :
		dimension(0), resolution(0), root(NULL)
	{
		root = new OctNode();
	}

	PtrOctTree(float dim, float res) :
		dimension(dim), resolution(res)
	{
		root = new OctNode();
	}

	PtrOctTree(float dim, float res, const vector<MolSphere>& mol) :
		dimension(dim), resolution(res)
	{
		Cube cube(dimension);
		root = new OctNode();
		create(root, cube, mol);
	}

	PtrOctTree(const PtrOctTree& rhs) :
		dimension(rhs.dimension), resolution(rhs.resolution), root(NULL)
	{
		if (rhs.root)
		{
			root = new OctNode(*rhs.root);
		}
	}

	PtrOctTree& operator=(PtrOctTree rhs)
	{
		swap(*this, rhs); //rhs passed by value
		return *this;
	}

	virtual ~PtrOctTree()
	{
		if (root)
			delete root;
	}

	friend void swap(PtrOctTree& first, PtrOctTree& second)
	{
		// enable ADL (not necessary in our case, but good practice)
		using std::swap;

		swap(first.resolution, second.resolution);
		swap(first.dimension, second.dimension);
		swap(first.root, second.root);
	}

	//invert filled and unfilled
	void invert()
	{
		return root->invert();
	}

	//mogrifying intersection
	bool intersect(const PtrOctTree& rhs)
	{
		return root->intersect(rhs.root);
	}
	//mogrifying union
	bool unionWith(const PtrOctTree& rhs)
	{
		return root->unionWith(rhs.root);
	}

	//truncate (maximum included volume) to specified resolution
	//only nodes that are >= resolution and full are kept
	void truncate(double resolution)
	{
		root->truncate(resolution, dimension);
	}

	//expand (minimum surrounding volume) so that any non-empty nodes <= resolution
	//are made full
	void grow(double resolution)
	{
		root->grow(resolution, dimension);
	}

	//volume calculations that don't require creating a tmp tree
	float intersectVolume(const PtrOctTree& rhs) const
	{
		return root->intersectVolume(rhs.root, dimension);
	}
	float unionVolume(const PtrOctTree& rhs) const
	{
		return root->unionVolume(rhs.root, dimension);
	}

	//return total volume contained in octtree
	float volume() const
	{
		return root->volume(dimension);
	}

	//return number of leaves
	unsigned leaves() const
	{
		return root->leaves();
	}

	void clear()
	{
		root->deleteChildren();
	}
	void fill()
	{
		root->deleteChildren();
		root->type = Full;
	}

	void write(ostream& out) const;
	void read(istream& in);
};

#endif /* PTROCTTREE_H_ */
