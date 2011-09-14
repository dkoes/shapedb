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
	};

	struct GSSInternalNode;
	struct GSSNode
	{
		GSSInternalNode *parent;
		OctTree *MIV; //maximum included volume
		OctTree *MSV; //minimum surrounding volume
		float res; // resolution
		unsigned which; //which child this is

		GSSNode(float d, float r):parent(NULL), MIV(new OctTree(d,r)), MSV(new OctTree(d,r)), res(r),which(0) {

		}
		virtual ~GSSNode();

	};

	struct GSSInternalNode: public GSSNode
	{

		vector<GSSNode*> children;

		GSSInternalNode(float d, float r): GSSNode(d, r)
		{
			children.reserve(MaxSplit);
		}

		~GSSInternalNode()
		{
			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				delete children[i];
			}
			if(MIV) delete MIV;
			if(MSV) delete MSV;
		}

		void update(GSSTree& gTree, unsigned whichChild, GSSNode *newnode);

		void addChild(GSSNode *child);
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

		~GSSLeafNode()
		{
			for(unsigned i = 0, n = trees.size(); i < n; i++)
			{
				if(trees[i]) delete trees[i];
			}
		}

		void insert(GSSTree& gTree, const OctTree& tree, const LeafData& data);

	};

	float maxres;
	float min[3];
	float dim;

	void setBoundingBox(array<float,6>& box);

	GSSInternalNode *root;

	static const unsigned MaxSplit;

	static void split(const vector<OctTree*>& MIV, const vector<OctTree*>& MSV,
			vector<unsigned>& s1, vector<unsigned>& s2);
	void createRoot(GSSNode *left, GSSNode *right);



public:
	GSSTree(filesystem::path& in); //initialize into memory from file; may mmap

	//set the global min and max extents of the GSS tree and
	//the highest resolution (default 1A)
	GSSTree(array<float,6>& boundingbox, double maxr=1.0): maxres(maxr), root(NULL)
	{
		setBoundingBox(boundingbox);
	}

	virtual ~GSSTree() { delete root; root = NULL; }

	//add a single mol
	void add(const vector<MolSphere>& mol);

	//nearest neighbor search, return closest set of molspheres
	void nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res);

	void write(const filesystem::path& out); //dump to file
};

#endif /* GSSTREE_H_ */
