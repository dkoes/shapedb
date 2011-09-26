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

	struct GSSNode;
	struct GSSInternalNode;
	struct GSSLeafNode;

	class PartitionData
	{
	public:
		PartitionData();
		virtual unsigned size() const = 0;

		virtual const LinearOctTree* getMSV(unsigned i) const = 0;
		virtual const LinearOctTree* getMIV(unsigned i) const = 0;
	};

	class LeafPartitionData: public PartitionData
	{
		const vector<LinearOctTree*> *trees;
		const vector<LeafData> *data;
	public:

		LeafPartitionData(): trees(NULL), data(NULL) {}
		LeafPartitionData(const vector<LinearOctTree*> *Ts, const vector<LeafData> *Ds): trees(Ts), data(Ds)
		{
			assert(trees->size() == data->size());
		}

		LinearOctTree* getTree(unsigned i)
		{
			return (*trees)[i];
		}

		const LeafData& getData(unsigned i)
		{
			return (*data)[i];
		}

		virtual const LinearOctTree* getMSV(unsigned i) const
		{
			return (*trees)[i];
		}

		virtual const LinearOctTree* getMIV(unsigned i) const
		{
			return (*trees)[i];
		}

		unsigned size() const { return trees != NULL ? trees->size() : 0; }

	};

	class NodePartitionData: public LeafPartitionData
	{
		const vector<GSSNode*> *nodes;
	public:
		NodePartitionData(): nodes(NULL) {}
		NodePartitionData(const vector<GSSNode*> *Ns): nodes(Ns) {  }

		GSSNode* getNode(unsigned i) { return (*nodes)[i]; }

		virtual const LinearOctTree* getMSV(unsigned i) const
		{
			return (*nodes)[i]->MSV;
		}

		virtual const LinearOctTree* getMIV(unsigned i) const
		{
			return (*nodes)[i]->MIV;
		}
	};

	/*
	 * A  partioner is an interface into a collection of MIV/MSV trees
	 * that presents the same interface to the partitioner regardless of whether
	 * we're partitioning leaves or internal nodes
	 */
	class Partitioner
	{
	protected:
		vector<unsigned> tindex;
		const PartitionData *partdata;
		void init(unsigned n);
	public:

		Partitioner(): partdata(NULL) {}

		//initialize initial partition
		Partitioner(const PartitionData *part): partdata(part)
		{
			tindex.reserve(partdata->size());
			for(unsigned i = 0, n = partdata->size(); i < n; i++)
				tindex.push_back(i);
		}

		unsigned getDataIndex(unsigned i) const { return tindex[i]; }
		unsigned size() const { return tindex.size(); }
		bool findOctant(unsigned level, bool splitMSV, vector<unsigned>& octantcoord);
		void partitionOnOctant(const vector<unsigned>& octantcoord, bool splitMSV, vector<Partitioner>& parts);
		void add(const Partitioner& from);
		void addSingle(const Partitioner& from, unsigned fromindex);
		void packClusters(unsigned max, vector<Partitioner>& clusters);
		double shannonEntropy(unsigned *patternCnts, unsigned total);
		void getOctants(unsigned level, bool splitMSV, vector< vector<unsigned> >& octants);

	};


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
		LinearOctTree *MIV; //maximum included volume
		LinearOctTree *MSV; //minimum surrounding volume
		float res; // resolution
		unsigned which; //which child this is

		GSSNode(): parent(NULL), MIV(new LinearOctTree(0,0)), MSV(new LinearOctTree(0,0)), res(0), which(0)
		{
			MIV->fill();
		}
		GSSNode(float d, float r):parent(NULL), MIV(new LinearOctTree(d,r)), MSV(new LinearOctTree(d,r)), res(r),which(0) {
			MIV->fill();
		}
		virtual ~GSSNode();

		virtual void selfUpdate() = 0;
		void fullUpdate();
		//examine every leaf to find nearest - for debugging and testing
		virtual void scanNearest(const LinearOctTree& tree, float& distance, LeafData& data) = 0;
		//find "closest" object to tree and put data into data
		virtual void findNearest(const LinearOctTree& tree, float& distance, LeafData& data) = 0;

		virtual void scanTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res) = 0;
		virtual void findTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res) = 0;

		virtual void findInsertionPoint(const LinearOctTree& tree, float& distance, GSSLeafNode*& leaf) = 0;
		//find k-best leaves for tree
		virtual void findInsertionPoints(const LinearOctTree& tree, vector<LeafDistPair>& kbest, unsigned k) = 0;

		virtual void printLeafInfo(unsigned depth) const = 0;
		virtual unsigned size() const = 0;
		virtual unsigned numLeaves() const = 0;

		float combinedVolumeChange(LinearOctTree *miv, LinearOctTree *msv) const;

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
		void selfUpdate();
		virtual void scanNearest(const LinearOctTree& tree, float& distance, LeafData& data);
		virtual void findNearest(const LinearOctTree& tree, float& distance, LeafData& data);
		virtual void scanTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res);
		virtual void findTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res);
		virtual void findInsertionPoint(const LinearOctTree& tree, float& distance, GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const LinearOctTree& tree, vector<LeafDistPair>& kbest, unsigned k);

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
		vector<LinearOctTree*> trees;
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

		void selfUpdate();
		void insert(GSSTree& gTree, const LinearOctTree& tree, const LeafData& data, vector<LeafDistPair>& kbest);
		virtual void scanNearest(const LinearOctTree& tree, float& distance, LeafData& data);
		virtual void findNearest(const LinearOctTree& tree, float& distance, LeafData& data);
		virtual void scanTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res);
		virtual void findTweeners(const LinearOctTree& min, const LinearOctTree& max, vector<LeafData>& res);
		virtual void findInsertionPoint(const LinearOctTree& tree, float& distance, GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const LinearOctTree& tree, vector<LeafDistPair>& kbest, unsigned k);

		void moveTreeFrom(GSSLeafNode* from, unsigned t);

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
	unsigned maxlevel;
	GSSNode *root;

	static const unsigned MaxSplit;

	void setBoundingBox(array<float,6>& box);

	float deltaFit(const GSSLeafNode* to, GSSLeafNode* from, unsigned t);
	static bool fitsInbetween(const LinearOctTree *MIV, const LinearOctTree *MSV, const LinearOctTree *min, const LinearOctTree *max);
	static float leafDist(const LinearOctTree* obj, const LinearOctTree *leaf);
	static float searchDist(const LinearOctTree* obj, const LinearOctTree *MIV, const LinearOctTree *MSV, float& min, float& max);
	static float splitDist(const LinearOctTree* leftMIV, const LinearOctTree* leftMSV, const LinearOctTree* rightMIV, const LinearOctTree* rightMSV);
	static void split(const vector<LinearOctTree*>& MIV, const vector<LinearOctTree*>& MSV,
			vector<unsigned>& s1, vector<unsigned>& s2);
	void createRoot(GSSNode *left, GSSNode *right);

	void transformMol(const vector<MolSphere>& inspheres, vector<MolSphere>& outspheres);

	GSSLeafNode* leafFromPartition(Partitioner& partitioner, LeafPartitionData& leafdata);
	void partitionLeaves(Partitioner& partitioner, LeafPartitionData& leafdata, unsigned level, bool splitMSV, vector<GSSNode*>& nodes);
	GSSInternalNode* nodeFromPartition(Partitioner& partitioner, NodePartitionData& nodedata);
	void generateNextLevel(Partitioner& partitioner, NodePartitionData& nodedata,unsigned level, bool splitMSV,  vector<GSSNode*>& nodes);


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
	//bulk load a collection of molecules into a new tree
	void load(const vector<vector<MolSphere> >& mols);

	//nearest neighbor search, return closest set of molspheres
	void nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res);
	//distance constrained search - inbetween little and big
	void dc_search(const vector<MolSphere>& little, const vector<MolSphere>& big, vector<vector<MolSphere> >& res);
	//inclusion/exclusion search with explicit volumes
	void inex_search(const vector<MolSphere>& inc, const vector<MolSphere>& exc, vector<vector<MolSphere> >& res);
	//constrained search with trees as input
	void tree_range_search(const LinearOctTree& smallTree, const LinearOctTree& bigTree, vector<vector<MolSphere>  >& res);

	void write(ostream& out); //dump to file
	void read(istream& in); //dump to file

	void printRootInfo() const { cout << root->MIV->volume() <<","<<root->MSV->volume() << "\n"; } //for debug
	void printLeafInfo() const { root->printLeafInfo(0); cout << "\n"; } //for debugging
	unsigned size() const { return root->size(); }
	unsigned numLeaves() const { return root->numLeaves(); }
};

#endif /* GSSTREE_H_ */
