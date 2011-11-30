/*
 * GSSTreeSearcher.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *      Class for loading and searching a GSSTree.
 */

#ifndef GSSTREESEARCHER_H_
#define GSSTREESEARCHER_H_

#include "GSSTreeStructures.h"
#include "Molecule.h"
#include "MemMapped.h"
#include "MappableOctTree.h"


using namespace boost;
using namespace std;

typedef Molecule Object; //eventually template this


class GSSTreeSearcher
{
	MemMapped objects; //memory mapped objects
	MemMapped internalNodes;
	MemMapped leaves;

	bool verbose; //for debugging
	unsigned total;
	float dimension;
	float resolution;

	void findTweeners(const GSSInternalNode* node,  const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& res,unsigned level);
	void findTweeners(const GSSLeaf* node, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& res);

	struct ObjDist
	{
		file_index objpos;
		double dist;

		bool operator<(const ObjDist& rhs) const
		{
			return dist < rhs.dist;
		}
	};

	//maintain K best objects found
	class TopKObj
	{
		unsigned k;
		vector<ObjDist> objs;
	public:
		TopKObj(unsigned K): k(K)
		{
			objs.reserve(K+1);
		}

		void add(file_index pos, double dist);
		double worst() const
		{
			if(objs.size() < k) return HUGE_VAL;
			else return objs.back().dist;
		}

		unsigned size() const { return objs.size(); }
		const ObjDist& operator[](unsigned i) const { return objs[i]; }
	};

	void findNearest(const GSSInternalNode* node,  const MappableOctTree* obj, TopKObj& res,unsigned level);
	void findNearest(const GSSLeaf* node, const MappableOctTree* obj, TopKObj& res);

	unsigned fitsCheck;
	unsigned nodesVisited;
	vector<unsigned> levelCnts;
	vector<unsigned> usefulLevelCnts;
	vector<unsigned> maxlevelCnts;
	unsigned leavesVisited;
	unsigned fullLeaves;
	bool fitsInbetween(const MappableOctTree *MIV, const MappableOctTree *MSV,
			const MappableOctTree  *min, const MappableOctTree *max);
public:
	GSSTreeSearcher(bool v = false): verbose(v), total(0) {}

	bool load(const filesystem::path& dbpath);
	void clear();
	~GSSTreeSearcher();

	unsigned size() const { return total; }

	//return everything with a shape between smallObj and bigObj
	void dc_search(const Object& smallObj, const Object& bigObj, bool invertBig,
			vector<Object>& res);

	//linear scan
	void dc_scan_search(const Object& smallObj, const Object& bigObj, bool invertBig,
			vector<Object>& res);

	//return k objects closest to obj
	void nn_search(const Object& obj, unsigned k, vector<Object>& res);

	float getDimension() const { return dimension; }
	float getResolution() const { return resolution; }
};

#endif /* GSSTREESEARCHER_H_ */
