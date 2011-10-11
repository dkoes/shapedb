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
#include <climits>

#include "MolSphere.h"
#include "OctTree.h"
#include "OctTreeFactory.h"

using namespace boost;
using namespace std;

//the gss tree
class GSSTree
{
	struct LeafData
	{
		vector<MolSphere> spheres;

		LeafData()
		{
		}
		LeafData(const vector<MolSphere>& sph) :
			spheres(sph)
		{
		}

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

		virtual const OctTree* getMSV(unsigned i) const = 0;
		virtual const OctTree* getMIV(unsigned i) const = 0;
	};

	class LeafPartitionData: public PartitionData
	{
		const vector<OctTree*> *trees;
		const vector<LeafData> *data;
	public:

		LeafPartitionData() :
			trees(NULL), data(NULL)
		{
		}
		LeafPartitionData(const vector<OctTree*> *Ts,
				const vector<LeafData> *Ds) :
			trees(Ts), data(Ds)
		{
			assert(trees->size() == data->size());
		}

		OctTree* getTree(unsigned i)
		{
			return (*trees)[i];
		}

		const LeafData& getData(unsigned i)
		{
			return (*data)[i];
		}

		virtual const OctTree* getMSV(unsigned i) const
		{
			return (*trees)[i];
		}

		virtual const OctTree* getMIV(unsigned i) const
		{
			return (*trees)[i];
		}

		unsigned size() const
		{
			return trees != NULL ? trees->size() : 0;
		}

	};

	class NodePartitionData: public PartitionData
	{
		const vector<GSSNode*> *nodes;
	public:
		NodePartitionData() :
			nodes(NULL)
		{
		}
		NodePartitionData(const vector<GSSNode*> *Ns) :
			nodes(Ns)
		{
		}

		GSSNode* getNode(unsigned i)
		{
			return (*nodes)[i];
		}

		virtual const OctTree* getMSV(unsigned i) const
		{
			return (*nodes)[i]->MSV;
		}

		virtual const OctTree* getMIV(unsigned i) const
		{
			return (*nodes)[i]->MIV;
		}

		unsigned size() const
		{
			return nodes != NULL ? nodes->size() : 0;
		}
	};

	/*
	 * A  partioner is an interface into a collection of MIV/MSV trees
	 * that presents the same interface to the partitioner regardless of whether
	 * we're partitioning leaves or internal nodes
	 */
	class Partitioner
	{
	public:
		typedef vector<vector<unsigned> > OctantCoords;

	protected:
		vector<unsigned> tindex;
		const PartitionData *partdata;
		OctantCoords MSVoctants;
		OctantCoords MIVoctants;

		vector<bool> MSVdoneOct;
		vector<bool> MIVdoneOct;

		OctTree *MIV;
		OctTree *MSV;

		unsigned maxlevel;
		unsigned level;
		bool splitMSV;

		void updateOctantsForLevel(unsigned level, OctantCoords& octants, vector<bool>& done);
		void mergeWith(const Partitioner& rhs);

	public:

		Partitioner() :
			partdata(NULL), MIV(NULL), MSV(NULL), maxlevel(0), level(0), splitMSV(true)
		{
		}

		//initialize initial partition
		Partitioner(const PartitionData *part, unsigned maxl) :
			partdata(part), MIV(NULL), MSV(NULL), maxlevel(maxl), level(0), splitMSV(true)
		{

		}

		Partitioner(const Partitioner& rhs): tindex(rhs.tindex), partdata(rhs.partdata),
				MSVoctants(rhs.MSVoctants), MIVoctants(rhs.MIVoctants), MSVdoneOct(rhs.MSVdoneOct),
				MIVdoneOct(rhs.MIVdoneOct), MIV(NULL), MSV(NULL), maxlevel(rhs.maxlevel),level(rhs.level), splitMSV(rhs.splitMSV)
		{
			if(rhs.MIV)
				MIV = rhs.MIV->clone();
			if(rhs.MSV)
				MSV = rhs.MSV->clone();
		}

		Partitioner& operator=(Partitioner rhs)
		{
			swap(*this, rhs); //rhs passed by value
			return *this;
		}

		friend void swap(Partitioner& first, Partitioner& second)
		{
			// enable ADL (not necessary in our case, but good practice)
			using std::swap;

			swap(first.tindex,second.tindex);
			swap(first.partdata,second.partdata);
			swap(first.MSVoctants,second.MSVoctants);
			swap(first.MIVoctants,second.MIVoctants);

			swap(first.MSVdoneOct,second.MSVdoneOct);
			swap(first.MIVdoneOct,second.MIVdoneOct);

			swap(first.MIV,second.MIV);
			swap(first.MSV,second.MSV);
			swap(first.level, second.level);
			swap(first.splitMSV, second.splitMSV);

			swap(first.maxlevel,second.maxlevel);
		}


		~Partitioner()
		{
			if(MIV) delete MIV;
			if(MSV) delete MSV;
		}

		void clear()
		{
			*this = Partitioner();
		}

		//create a sub partition of parent without (as yet) the actual indices
		void inheritFrom(const Partitioner& parent)
		{
			partdata = parent.partdata;
			MSVoctants = parent.MSVoctants;
			MIVoctants = parent.MIVoctants;
			MSVdoneOct = parent.MSVdoneOct;
			MIVdoneOct = parent.MIVdoneOct;
			tindex.reserve(parent.tindex.size());

			level = parent.level;
			splitMSV = parent.splitMSV;
			maxlevel = parent.maxlevel;
		}

		void initFromData()
		{
			tindex.reserve(partdata->size());
			for (unsigned i = 0, n = partdata->size(); i < n; i++)
			{
				tindex.push_back(i);
				if(MSV == NULL)
					MSV = partdata->getMSV(i)->clone();
				else
					MSV->unionWith(partdata->getMSV(i));

				if(MIV == NULL)
					MIV = partdata->getMIV(i)->clone();
				else
					MIV->intersect(partdata->getMIV(i));
			}
		}

		unsigned getDataIndex(unsigned i) const
		{
			return tindex[i];
		}
		unsigned size() const
		{
			return tindex.size();
		}
		bool findOctantAndSetUsed(unsigned level, bool splitMSV, vector<unsigned>& octantcoord);
		void partitionOnOctant(const vector<unsigned>& octantcoord,
				bool splitMSV, vector<Partitioner>& parts);
		void add(const Partitioner& from);
		void addSingle(const Partitioner& from, unsigned fromindex);
		void packClusters(unsigned max, vector<Partitioner>& clusters);
		void kClusters(unsigned k, vector<Partitioner>& clusters);
		double shannonEntropy(unsigned *patternCnts, unsigned total);

		const OctTree* getMSV() const { return MSV; }
		const OctTree* getMIV() const { return MIV; }

		static void clusterPartitions(vector<Partitioner>& clusters);

		void topDownOctantPartition(vector<Partitioner>& parts);
		void topDownKSamplePartition(vector<Partitioner>& parts);

		bool unableToPartition() const;

		void getCenter(const OctTree *& MIV, const OctTree *& MSV) const;

	};

	//class for storing leaves with a distance so we can sort
	struct LeafDistPair
	{
		GSSLeafNode *leaf;
		float distance;

		LeafDistPair() :
			leaf(NULL), distance(HUGE_VAL)
		{
		}
		LeafDistPair(GSSLeafNode *l, float d) :
			leaf(l), distance(d)
		{
		}

		bool operator<(const LeafDistPair& rhs) const
		{
			return distance < rhs.distance;
		}
	};
	struct Stats
	{
		unsigned min;
		unsigned max;
		unsigned total;
		unsigned cnt;
		unsigned singletonCnt;

		Stats() :
			min(UINT_MAX), max(0), total(0), cnt(0), singletonCnt(0)
		{
		}
	};

	struct GSSNode
	{
		GSSInternalNode *parent;
		OctTree *MIV; //maximum included volume
		OctTree *MSV; //minimum surrounding volume
		float res; // resolution
		unsigned which; //which child this is

		GSSNode(const OctTreeFactory& octGen) :
			parent(NULL), MIV(octGen.newOctTree(0, 0)),
					MSV(octGen.newOctTree(0, 0)), res(0), which(0)
		{
			MIV->fill();
		}
		GSSNode(const OctTreeFactory& octGen, float d, float r) :
			parent(NULL), MIV(octGen.newOctTree(d, r)),
					MSV(octGen.newOctTree(d, r)), res(r), which(0)
		{
			MIV->fill();
		}
		virtual ~GSSNode();

		virtual void selfUpdate() = 0;
		void fullUpdate();
		//examine every leaf to find nearest - for debugging and testing
		virtual void scanNearest(const OctTree* tree, float& distance,
				LeafData& data) = 0;
		//find "closest" object to tree and put data into data
		virtual void findNearest(const OctTree* tree, float& distance,
				LeafData& data) = 0;

		virtual void scanTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res) = 0;
		virtual void findTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res) = 0;

		virtual void findInsertionPoint(const OctTree*tree, float& distance,
				GSSLeafNode*& leaf) = 0;
		//find k-best leaves for tree
		virtual void findInsertionPoints(const OctTree* tree,
				vector<LeafDistPair>& kbest, unsigned k) = 0;

		virtual void printLeafInfo(unsigned depth) const = 0;
		virtual unsigned size() const = 0;
		virtual unsigned numLeaves() const = 0;

		float combinedVolumeChange(const OctTree *miv, const OctTree *msv) const;

		virtual void write(ostream& out) const;
		virtual void read(istream& in, GSSInternalNode *parent);
		static GSSNode* readCreate(const OctTreeFactory& octGen, istream& in,
				GSSInternalNode *parent);

		virtual unsigned getStats(Stats& leaves, Stats& nodes) const = 0;
	};

	struct GSSInternalNode: public GSSNode
	{

		vector<GSSNode*> children;

		GSSInternalNode(const OctTreeFactory& octGen, float d, float r) :
			GSSNode(octGen, d, r)
		{
			children.reserve(MaxSplit);
		}

		GSSInternalNode(const OctTreeFactory& octGen) :
			GSSNode(octGen)
		{
		}

		virtual ~GSSInternalNode()
		{
			for (unsigned i = 0, n = children.size(); i < n; i++)
			{
				delete children[i];
			}
		}

		void update(GSSTree& gTree, unsigned whichChild, GSSNode *newnode);
		void selfUpdate();
		virtual void scanNearest(const OctTree* tree, float& distance,
				LeafData& data);
		virtual void findNearest(const OctTree* tree, float& distance,
				LeafData& data);
		virtual void scanTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res);
		virtual void findTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res);
		virtual void findInsertionPoint(const OctTree* tree, float& distance,
				GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const OctTree* tree,
				vector<LeafDistPair>& kbest, unsigned k);

		void addChild(GSSNode *child);

		void printLeafInfo(unsigned depth) const
		{
			for (unsigned i = 0, n = children.size(); i < n; i++)
			{
				children[i]->printLeafInfo(depth + 1);
			}
		}

		unsigned size() const
		{
			unsigned ret = 0;
			for (unsigned i = 0, n = children.size(); i < n; i++)
			{
				ret += children[i]->size();
			}
			return ret;
		}

		unsigned numLeaves() const
		{
			unsigned ret = 0;
			for (unsigned i = 0, n = children.size(); i < n; i++)
			{
				ret += children[i]->numLeaves();
			}
			return ret;
		}

		virtual void write(ostream& out) const;
		virtual void read(const OctTreeFactory& octGen, istream& in,
				GSSInternalNode *parPtr);

		unsigned getStats(Stats& leaves, Stats& nodes) const;

	};

	struct GSSLeafNode: public GSSNode
	{
		vector<OctTree*> trees;
		vector<LeafData> data;

		GSSLeafNode(const OctTreeFactory& octGen, float d, float r) :
			GSSNode(octGen, d, r)
		{
			trees.reserve(MaxSplit);
			data.reserve(MaxSplit);
		}

		GSSLeafNode(const OctTreeFactory& octGen) :
			GSSNode(octGen)
		{
		}

		~GSSLeafNode()
		{
			for (unsigned i = 0, n = trees.size(); i < n; i++)
			{
				if (trees[i])
					delete trees[i];
			}
		}

		void selfUpdate();
		void insert(GSSTree& gTree, OctTree* tree, const LeafData& data,
				vector<LeafDistPair>& kbest);
		virtual void scanNearest(const OctTree* tree, float& distance,
				LeafData& data);
		virtual void findNearest(const OctTree* tree, float& distance,
				LeafData& data);
		virtual void scanTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res);
		virtual void findTweeners(const OctTree* min, const OctTree* max,
				vector<LeafData>& res);
		virtual void findInsertionPoint(const OctTree* tree, float& distance,
				GSSLeafNode*& leaf);
		virtual void findInsertionPoints(const OctTree* tree,
				vector<LeafDistPair>& kbest, unsigned k);

		void moveTreeFrom(GSSLeafNode* from, unsigned t);

		void printLeafInfo(unsigned depth) const
		{
			cout << depth << "(" << MIV->volume() << "," << MSV->volume()
					<< "):" << trees.size() << " ";
		}

		unsigned size() const
		{
			return trees.size();
		}

		unsigned numLeaves() const
		{
			return 1;
		}

		virtual void write(ostream& out) const;
		virtual void read(const OctTreeFactory& octGen, istream& in,
				GSSInternalNode *parPtr);

		unsigned getStats(Stats& leaves, Stats& nodes) const;

	};

	float maxres;
	float min[3];
	float dim;
	unsigned maxlevel;
	GSSNode *root;

	OctTreeFactory octGen;
	static const unsigned MaxSplit;

	void setBoundingBox(array<float, 6>& box);

	float deltaFit(const GSSLeafNode* to, GSSLeafNode* from, unsigned t);
	static bool fitsInbetween(const OctTree *MIV, const OctTree *MSV,
			const OctTree *min, const OctTree *max);
	static float leafDist(const OctTree* obj, const OctTree *leaf);
	static float searchDist(const OctTree* obj, const OctTree *MIV,
			const OctTree *MSV, float& min, float& max);
	static float splitDist(const OctTree* leftMIV, const OctTree* leftMSV,
			const OctTree* rightMIV, const OctTree* rightMSV);
	static void split(const vector<OctTree*>& MIV, const vector<OctTree*>& MSV,
			vector<unsigned>& s1, vector<unsigned>& s2);
	void createRoot(GSSNode *left, GSSNode *right);

	void transformMol(const vector<MolSphere>& inspheres,
			vector<MolSphere>& outspheres);

	GSSLeafNode* leafFromPartition(Partitioner& partitioner,
			LeafPartitionData& leafdata);
	void partitionLeaves(Partitioner& partitioner, LeafPartitionData& leafdata,
			vector<GSSNode*>& nodes);
	GSSInternalNode* nodeFromPartition(Partitioner& partitioner,
			NodePartitionData& nodedata);
	void partitionNodes(Partitioner& partitioner, NodePartitionData& nodedata,
			vector<GSSNode*>& nodes);

public:

	//set the global min and max extents of the GSS tree and
	//the highest resolution (default 1A)
	GSSTree(const OctTreeFactory& octfact, array<float, 6>& boundingbox,
			double maxr = 1.0) :
		maxres(maxr), root(NULL), octGen(octfact)
	{
		setBoundingBox(boundingbox);
		root = new GSSLeafNode(octGen, dim, maxres);
	}

	GSSTree() :
		maxres(0), root(NULL)
	{
		//basically can just read in
	}

	virtual ~GSSTree()
	{
		delete root;
		root = NULL;
	}

	//add a single mol
	void add(const vector<MolSphere>& mol);
	//bulk load a collection of molecules into a new tree
	void load(const vector<vector<MolSphere> >& mols);

	//nearest neighbor search, return closest set of molspheres
	void nn_search(const vector<MolSphere>& mol, vector<MolSphere>& res);
	//distance constrained search - inbetween little and big
	void dc_search(const vector<MolSphere>& little,
			const vector<MolSphere>& big, vector<vector<MolSphere> >& res);
	//inclusion/exclusion search with explicit volumes
	void inex_search(const vector<MolSphere>& inc,
			const vector<MolSphere>& exc, vector<vector<MolSphere> >& res);
	//constrained search with trees as input
	void tree_range_search(const OctTree* smallTree, const OctTree* bigTree,
			vector<vector<MolSphere> >& res);

	void write(ostream& out); //dump to file
	void read(istream& in); //dump to file

	void printRootInfo() const
	{
		cout << root->MIV->volume() << "," << root->MSV->volume() << "\n";
	} //for debug
	void printLeafInfo() const
	{
		root->printLeafInfo(0);
		cout << "\n";
	} //for debugging
	unsigned size() const
	{
		return root->size();
	}
	unsigned numLeaves() const
	{
		return root->numLeaves();
	}

	void printStats() const; //for debugging
};

#endif /* GSSTREE_H_ */
