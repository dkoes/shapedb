/*
 * GSSTree.h
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 *
 *      http://dx.doi.org/10.1145/304181.304219
 *
 *      A GSS tree for molecular data.  When building takes as input
 *      a vector of spheres (x,y,z,r) that represent a molecule
 *      (eventaully will add meta data).
 */

#ifndef GSSTREE_H_
#define GSSTREE_H_

#include <boost/filesystem.hpp>
#include <boost/array.hpp>
#include <vector>
#include <iostream>
#include <cmath>

#include "MolSphere.h"
#include "OctTree.h"

using namespace boost;
using namespace std;



//the gss tree
class GSSTree
{
	struct LeafData
	{
		vector<MolSphere> spheres;

		LeafData() {}
		LeafData(const vector<MolSphere>& sph): spheres(sph) {}

		void write(ostream& out) const;
		void read(istream& in);
	};

	struct GSSInternalNode;
	struct GSSLeafNode;

	//class for storing leaves with a distance so we can sort
	struct LeafDistPair
	{
		GSSLeafNode *leaf;
		float distance;

		LeafDistPair(): leaf(NULL), distance(HUGE_VAL) {}
		LeafDistPair(GSSLeafNode *l, float d): leaf(l), distance(d) {}

		bool operator<(const LeafDistPair& rhs) const
		{
			return distance < rhs.distance;
		}
	};

	struct GSSNode
	{
		GSSInternalNode *parent;
		OctTree *MIV; //maximum included volume
		OctTree *MSV; //minimum surrounding volume
		float res; // resolution
		unsigned which; //which child this is

		GSSNode(): parent(NULL), MIV(new OctTree(0,0)), MSV(new OctTree(0,0)), res(0), which(0)
		{
			MIV->fill();
		}
		GSSNode(float d, float r):parent(NULL), MIV(new OctTree(d,r)), MSV(new OctTree(d,r)), res(r),which(0) {
			MIV->fill();
		}
		virtual ~GSSNode();

		//examine every leaf to find nearest - for debugging and testing
		virtual void scanNearest(const OctTree& tree, float& distance, LeafData& data) = 0;
		//find "closest" object to tree and put data into data
		virtual void findNearest(const OctTree& tree, float& distance, LeafData& data) = 0;

		virtual void findInsertionPoint(const OctTree& tree, float& distance, GSSLeafNode*& leaf) = 0;
		//find k-best leaves for tree
		virtual void findInsertionPoints(const OctTree& tree, vector<LeafDistPair>& kbest, unsigned k) = 0;

		virtual void printLeafInfo(unsigned depth) const = 0;
		virtual unsigned size() const = 0;
		virtual unsigned numLeaves() const = 0;

		float combinedVolumeChange(OctTree *miv, OctTree *msv) const;

		virtual void write(ostream& out) const;
		virtual void read(istream& in, GSSInternalNode *parent);
		static GSSNode* readCreate(istream& in, GSSInternalNode *parent);

	};

	struct GSSInternalNode: public GSSNode
	{

		vector<GSSNode*> children;

		GSSInternalNode(float d, float r): GSSNode(d, r)
		{
			children.reserve(MaxSplit);
		}

		GSSInternalNode() {}

		virtual ~GSSInternalNode()
		{
			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				delete children[i];
			}
		}

		void update(GSSTree& gTree, unsigned whichChild, GSSNode *newnode);
		virtual void scanNearest(const OctTree& tree, float& distance, LeafData& data);
		virtual void findNearest(const OctTree& tree, float& distance, LeafData& data);
		virtual void findInsertionPoint(const OctTree& tree, float& distance, GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const OctTree& tree, vector<LeafDistPair>& kbest, unsigned k);

		void addChild(GSSNode *child);

		void printLeafInfo(unsigned depth) const
		{
			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				children[i]->printLeafInfo(depth+1);
			}
		}

		unsigned size() const
		{
			unsigned ret = 0;
			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				ret += children[i]->size();
			}
			return ret;
		}

		unsigned numLeaves() const
		{
			unsigned ret = 0;
			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				ret += children[i]->numLeaves();
			}
			return ret;
		}

		virtual void write(ostream& out) const;
		virtual void read(istream& in, GSSInternalNode *parPtr);
	};

	struct GSSLeafNode: public GSSNode
	{
		vector<OctTree*> trees;
		vector<LeafData> data;

		GSSLeafNode(float d, float r): GSSNode(d, r)
		{
			trees.reserve(MaxSplit);
			data.reserve(MaxSplit);
		}

		GSSLeafNode() {}

		~GSSLeafNode()
		{
			for(unsigned i = 0, n = trees.size(); i < n; i++)
			{
				if(trees[i]) delete trees[i];
			}
		}

		void insert(GSSTree& gTree, const OctTree& tree, const LeafData& data);
		virtual void scanNearest(const OctTree& tree, float& distance, LeafData& data);
		virtual void findNearest(const OctTree& tree, float& distance, LeafData& data);

		virtual void findInsertionPoint(const OctTree& tree, float& distance, GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const OctTree& tree, vector<LeafDistPair>& kbest, unsigned k);

		void printLeafInfo(unsigned depth) const
		{
			cout << depth << "(" << MIV->volume() << "," << MSV->volume() << "):" << trees.size() << " ";
		}

		unsigned size() const
		{
			return trees.size();
		}

		unsigned numLeaves() const { return 1; }

		virtual void write(ostream& out) const;
		virtual void read(istream& in, GSSInternalNode *parPtr);

	};

	float maxres;
	float min[3];
	float dim;
	GSSNode *root;

	static const unsigned MaxSplit;

	void setBoundingBox(array<float,6>& box);
	static float leafDist(const OctTree* obj, const OctTree *leaf);
	static float searchDist(const OctTree* obj, const OctTree *MIV, const OctTree *MSV, float& min, float& max);
	static float splitDist(OctTree* leftMIV, OctTree* leftMSV, OctTree* rightMIV, OctTree* rightMSV);
	static void split(const vector<OctTree*>& MIV, const vector<OctTree*>& MSV,
			vector<unsigned>& s1, vector<unsigned>& s2);
	void createRoot(GSSNode *left, GSSNode *right);

	void transformMol(const vector<MolSphere>& inspheres, vector<MolSphere>& outspheres);

public:

	//set the global min and max extents of the GSS tree and
	//the highest resolution (default 1A)
	GSSTree(array<float,6>& boundingbox, double maxr=1.0): maxres(maxr), root(NULL)
	{
		setBoundingBox(boundingbox);
		root = new GSSLeafNode(dim, maxres);
	}

	GSSTree(): maxres(0), root(NULL)
	{
		//basically can just read in
	}

	virtual ~GSSTree() { delete root; root = NULL; }

	//add a single mol
	void add(const vector<MolSphere>& mol);

	//nearest neighbor search, return closest set of molspheres
	void nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res);

	void write(ostream& out); //dump to file
	void read(istream& in); //dump to file

	void printRootInfo() const { cout << root->MIV->volume() <<","<<root->MSV->volume() << "\n"; } //for debug
	void printLeafInfo() const { root->printLeafInfo(0); cout << "\n"; } //for debugging
	unsigned size() const { return root->size(); }
	unsigned numLeaves() const { return root->numLeaves(); }
};

#endif /* GSSTREE_H_ */
