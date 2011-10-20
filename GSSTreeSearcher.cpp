/*
 * GSSTreeSearcher.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#include "GSSTreeSearcher.h"
#include <sys/types.h>
#include <sys/mman.h>
#include "MappableOctTree.h"
#include "Timer.h"

//load a gsstree database by mmapping files, return true if successfull
bool GSSTreeSearcher::load(const filesystem::path& dbpath)
{
	clear();
	//read in info
	filesystem::path infofile = dbpath / "info";
	ifstream info(infofile.string().c_str());
	if(!info)
		return false;
	unsigned levels;

	info >> dimension  >> resolution >> levels  >> total;

	//memory map files form db directory
	nodes.resize(levels);
	for(unsigned i = 0; i < levels; i++)
	{
		stringstream name;
		name << "level" << i;
		filesystem::path level = dbpath/name.str();
		if(!nodes[i].map(level.string(), true, true))
		{
			clear();
			return false;
		}
	}

	filesystem::path objfile = dbpath / "objs";
	if(!objects.map(objfile.string(), true, false))
	{
		clear();
		return false;
	}

	return true;
}



void GSSTreeSearcher::clear()
{
	total = 0;
	objects.clear();
	for(unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		nodes[i].clear();
	}
	nodes.clear();
}

GSSTreeSearcher::~GSSTreeSearcher()
{
	clear();
}

//find all the shapes in the database that lie bewtween smallObj and bigObj
//if invertBig is set, than treat as an excluded volume
void GSSTreeSearcher::dc_search(const Object& smallObj, const Object& bigObj, bool invertBig,
		vector<Object>& res)
{
	res.clear();
	const MappableOctTree *smallTree = MappableOctTree::create(dimension, resolution, smallObj);
	MappableOctTree *bigTree = MappableOctTree::create(dimension, resolution, bigObj);

	if(invertBig)
		bigTree->invert(dimension);

	const GSSNodeCommon* root = (GSSNodeCommon*)nodes.back().begin();
	Timer t;
	vector<file_index> respos;
	findTweeners(root, nodes.size()-1, smallTree, bigTree, respos);

	//extract objects in sequential order
	sort(respos.begin(), respos.end());
	for(unsigned i = 0, n = respos.size(); i < n; i++)
	{
		const char * addr = objects.begin()+respos[i];
		res.push_back(Object(addr));
	}

	if(verbose)
	{
		cout << "Found " << res.size() << " objects out of " << total << " in " << t.elapsed() << "s\n";
	}
	delete smallTree;
	delete bigTree;
}


//return true if the object(s) represented by MIV/MSV might fit in between min and max
static bool fitsInbetween(const MappableOctTree *MIV, const MappableOctTree *MSV,
		const MappableOctTree  *min, const MappableOctTree *max)
{
	//the MSV must completely enclose min
	if (!min->containedIn(MSV))
		return false;
	//MIV must be completely enclosed by max
	if (!MIV->containedIn(max))
		return false;

	return true;
}


//find everyting between small and big using linear scan
void GSSTreeSearcher::dc_scan_search(const Object& smallObj, const Object& bigObj, bool invertBig,
		vector<Object>& res)
{
	res.clear();
	const MappableOctTree *smallTree = MappableOctTree::create(dimension, resolution, smallObj);
	MappableOctTree *bigTree = MappableOctTree::create(dimension, resolution, bigObj);
	if(invertBig)
		bigTree->invert(dimension);
	const GSSLeaf* leaf = (GSSLeaf*)nodes.front().begin();
	const GSSLeaf* end = (GSSLeaf*)nodes.front().end();

	Timer t;
	vector<file_index> respos;
	for( ; leaf != end; leaf = (const GSSLeaf*)((char*)leaf + leaf->bytes()))
	{
		findTweeners(leaf, smallTree, bigTree, respos);
	}

	//extract objects in sequential order
	sort(respos.begin(), respos.end());
	for(unsigned i = 0, n = respos.size(); i < n; i++)
	{
		const char * addr = objects.begin()+respos[i];
		res.push_back(Object(addr));
	}

	if(verbose)
	{
		cout << "Scanned " << res.size() << " objects out of " << total << " in " << t.elapsed() << "s\n";
	}
	delete smallTree;
	delete bigTree;
}


//dispatch between leaf and internal node
void GSSTreeSearcher::findTweeners(const GSSNodeCommon* node, unsigned level, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& respos)
{
	if(node->isLeaf)
		findTweeners((const GSSLeaf*)node, min, max, respos);
	else
		findTweeners((const GSSInternalNode*)node, level, min, max, respos);
}

void GSSTreeSearcher::findTweeners(const GSSInternalNode* node, unsigned level, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& respos)
{
	unsigned nextlevel = level-1;
	assert(level > 0); //otherwise shoudl be leaf
	for(unsigned i = 0, n = node->size(); i < n; i++)
	{
		const GSSInternalNode::Child *child = node->getChild(i);

		if(fitsInbetween(child->getMIV(), child->getMSV(), min, max))
		{
			const GSSNodeCommon* next = (const GSSNodeCommon*)(nodes[nextlevel].begin() + child->node_pos);
			findTweeners(next, nextlevel, min, max, respos);
		}
	}
}

//identify and trees in this leaf that fit
void GSSTreeSearcher::findTweeners(const GSSLeaf* node, const MappableOctTree* min, const MappableOctTree* max, vector<file_index>& respos)
{
	for(unsigned i = 0, n = node->size(); i < n; i++)
	{
		const GSSLeaf::Child *child = node->getChild(i);

		if(fitsInbetween(&child->tree, &child->tree, min, max))
		{
			respos.push_back(child->object_pos);
		}
	}
}


