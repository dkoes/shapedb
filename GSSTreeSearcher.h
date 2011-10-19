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
	vector<MemMapped> nodes;

	bool verbose; //for debugging
	unsigned total;
	float dimension;
	float resolution;

	void findTweeners(const GSSNodeCommon* node, unsigned level, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& res);
	void findTweeners(const GSSInternalNode* node, unsigned level, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& res);
	void findTweeners(const GSSLeaf* node, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& res);

public:
	GSSTreeSearcher(bool v = false): verbose(v), total(0) {}

	bool load(const filesystem::path& dbpath);
	void clear();
	~GSSTreeSearcher();

	unsigned size() const { return total; }

	//return everything with a shape between smallObj and bigObj
	void dc_search(const Object& smallObj, const Object& bigObj,
			vector<Object>& res);

	//linear scan
	void dc_scan_search(const Object& smallObj, const Object& bigObj,
			vector<Object>& res);

};

#endif /* GSSTREESEARCHER_H_ */
