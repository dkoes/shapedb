/*
 * LinearOctTree.h
 *
 * This is an implementation of an oct tree that just stores a list
 * of the tree leaves with the interior nodes implicit in the arrangement
 * of the leaves.  Very space efficient, but indexing down to a specific
 * octant is O(n)
 *
 *  Created on: Sep 28, 2011
 *      Author: dkoes
 */

#ifndef LINEAROCTTREE_H_
#define LINEAROCTTREE_H_

#include "OctTree.h"

/* linearlized oct tree - tree is represented as a vector of values
 * we manage to not store children nodes this way, but always have to
 * iterate over the whole thing (no shortcuts)
 * */
class LinearOctTree: public OctTree
{
	float dimension;
	float resolution;

	enum Type
	{
		Empty, Full
	};
	struct OctVal
	{
		Type flag :1;
		unsigned level :7;

		OctVal() :
			flag(Empty), level(0)
		{
		}
		OctVal(Type t, unsigned l) :
			flag(t), level(l)
		{
		}
	}__attribute__((__packed__));

	vector<OctVal> tree;
	//store the volume of a cube at a given level
	vector<float> levelVolumes;

	void setLevelVolumes();
	//generate a linear oct tree from a mol
	void create(const Cube& cube, unsigned level, const vector<MolSphere>& mol);

	//shared code of intersection and union
	bool operation(const LinearOctTree& rhs, bool doUnion);
	float volOperation(const LinearOctTree& rhs, bool doUnion) const;

	static unsigned absorbTreeAtLevel(const vector<OctVal>& T, unsigned pos,
			unsigned level);
	static unsigned appendTreeAtLevel(const vector<OctVal>& T, unsigned pos,
			unsigned level, vector<OctVal>& appendto);
	unsigned volumeOfTreeAtLevel(const vector<OctVal>& T, unsigned pos,
			unsigned level, float& vol) const;
	static unsigned traverseOctantCoord(const vector<OctVal>& T,
			const vector<unsigned>& coord);

	//class for iterating over leaves of a tree, maining the cube of the leaf
	class LeafCubeIterator
	{
		unsigned pos;
		unsigned curlevel;
		const vector<OctVal>& tree;

		vector<Cube> nestedCubes;
		vector<unsigned> nestedOctants;
		void step(); //go to next leaf
	public:
		LeafCubeIterator(const LinearOctTree& oct):  pos(0), curlevel(0), tree(oct.tree)
		{
			nestedCubes.push_back(Cube(oct.dimension));
			while(curlevel < tree[0].level)
			{
				nestedCubes.push_back(nestedCubes.back().getOctant(0));
				nestedOctants.push_back(0);
				curlevel++;
			}

			while(tree[pos].flag == Empty)
				step();
		}

		operator bool() const //iterator valid
		{
			return pos < tree.size();
		}

		LeafCubeIterator& operator++()
		{
			step();
			return *this;
		}

		Cube operator*() const
		{
			return nestedCubes.back();
		}

	};

public:
	LinearOctTree() :
		dimension(0), resolution(0)
	{
		clear();
	}

	LinearOctTree(float dim, float res) :
		dimension(dim), resolution(res)
	{
		setLevelVolumes();
		clear();
	}

	LinearOctTree(float dim, float res, const vector<MolSphere>& mol) :
		dimension(dim), resolution(res)
	{
		setLevelVolumes();
		Cube cube(dimension);
		create(cube, 0, mol);
	}

	LinearOctTree(const LinearOctTree& rhs) :
		dimension(rhs.dimension), resolution(rhs.resolution), tree(rhs.tree),
				levelVolumes(rhs.levelVolumes)
	{
	}

	virtual OctTree* clone() const
	{
		return new LinearOctTree(*this);
	}

	LinearOctTree& operator=(LinearOctTree rhs)
	{
		swap(*this, rhs); //rhs passed by value
		return *this;
	}

	virtual ~LinearOctTree()
	{
	}

	friend void swap(LinearOctTree& first, LinearOctTree& second)
	{
		// enable ADL (not necessary in our case, but good practice)
		using std::swap;

		swap(first.resolution, second.resolution);
		swap(first.dimension, second.dimension);
		swap(first.tree, second.tree);
		swap(first.levelVolumes, second.levelVolumes);
	}

	//invert filled and unfilled
	void invert();

	//mogrifying intersection
	bool intersect(const OctTree* rhs);
	//mogrifying union
	bool unionWith(const OctTree* rhs);

	//volume calculations that don't require creating a tmp tree
	float intersectVolume(const OctTree* rhs) const;
	float unionVolume(const OctTree* rhs) const;

	//return total volume contained in octtree
	float volume() const;

	//return number of leaves
	unsigned leaves() const
	{
		return tree.size();
	}

	void clear()
	{
		tree.clear();
		tree.push_back(OctVal(Empty, 0));
	}
	void fill()
	{
		tree.clear();
		tree.push_back(OctVal(Full, 0));
	}

	void write(ostream& out) const;
	void read(istream& in);

	unsigned getOctantPattern(const vector<unsigned>& coord, bool MSV) const;

	unsigned countInteriorNodesAtLevel(unsigned level) const;

	float hausdorffDistance(const OctTree* B) const;
};

#endif /* LINEAROCTTREE_H_ */
