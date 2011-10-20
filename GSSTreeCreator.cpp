/*
 * GSSTreeCreator.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "GSSTreeCreator.h"
#include "DataViewers.h"
#include "TopDownPartitioner.h"


//convience function for creating an indexed path name
static string nextString(filesystem::path p, const char *base, unsigned i)
{
	stringstream str;
	str << base << i;
	return filesystem::path(p/ str.str()).string();
}

//return true if successfull
bool GSSTreeCreator::create(filesystem::path dir, Object::iterator& itr, float dim, float res)
{
	WorkFile currenttrees;
	WorkFile nexttrees;
	//create directory
	if(filesystem::exists(dir))
	{
		cerr << dir << " already exists.  Exiting\n";
		return false;
	}
	if (!filesystem::create_directory(dir))
	{
		cerr << "Unable to create database directory ";
		return false;
	}
	dbpath = dir;

	filesystem::path objfile = dbpath / "objs";
	string curtreesfile = filesystem::path(dbpath / "trees").string();
	string nexttreesfile = filesystem::path(dbpath / "nexttrees").string();

	//write out objects and trees
	objects.set(objfile.string().c_str());
	currenttrees.set(curtreesfile.c_str());
	vector<file_index> treeindices;
	vector<file_index> objindices;
	unsigned cnt = 0;
	for( ; itr ; ++itr)
	{
		const Object& obj = *itr;
		objindices.push_back((file_index)objects.file->tellp());
		obj.write(*objects.file);

		//leaf object
		treeindices.push_back((file_index)currenttrees.file->tellp());
		MappableOctTree *tree = MappableOctTree::create(dim, res, obj);
		tree->write(*currenttrees.file);
		delete tree;
		cnt++;
	}

	//partition leaves into bottom level
	//setup level
	nodes.push_back(WorkFile(nextString(dbpath, "level", nodes.size()).c_str()));

	vector<file_index> nodeindices; nodeindices.reserve(objindices.size()/2);

	//setup next trees
	nexttrees.set(nexttreesfile.c_str());

	//map the data
	currenttrees.switchToMap();
	//this clears the indices
	LeafViewer leafdata(currenttrees.map->get_address(), treeindices, objindices);
	leveler->createNextLevel(leafdata, nodes.back().file, nodeindices, nexttrees.file, treeindices);

	while(nodeindices.size() > 1)
	{
		//setup curr/next trees
		currenttrees.remove();
		swap(currenttrees, nexttrees);
		currenttrees.switchToMap();

		//this stores and resets the indices
		NodeViewer nodedata(currenttrees.map->get_address(), treeindices, nodeindices);

		//setup next level
		nodes.push_back(WorkFile(nextString(dbpath, "level", nodes.size()).c_str()));

		nexttrees.set(nextString(dbpath, "trees", nodes.size()).c_str());

		leveler->createNextLevel(nodedata, nodes.back().file, nodeindices, nexttrees.file, treeindices);
	}
	currenttrees.remove();
	nexttrees.remove();

	optimizeLevels();

	//output general info
	filesystem::path infoname = dbpath / "info";
	ofstream info(infoname.string().c_str());
	info << dim <<" " << res <<" " << nodes.size() << " " << cnt << "\n";

	//clear workfile memory
	for(unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		nodes[i].clear();
	}
	return true;
}


//recursive helper for optimizing level output
void GSSTreeCreator::optimizeLevelsR(ostream& outnodes, ostream& outleaves, const GSSNodeCommon *n, unsigned level)
{
	if(n->isLeaf)
	{
		const GSSLeaf* leaf = (const GSSLeaf*)n;
		outleaves.write((const char*)leaf, leaf->bytes());
	}
	else
	{
		const GSSInternalNode* node = (const GSSInternalNode*)n;
		outnodes.write((const char*)node, node->bytes());
		unsigned nextlevel = level-1;

		for(unsigned i = 0, n = node->size(); i < n; i++)
		{
			const GSSInternalNode::Child *child = node->getChild(i);
			const GSSNodeCommon* next = (const GSSNodeCommon*)((const char*)nodes[nextlevel].map->get_address() + child->node_pos);
			optimizeLevelsR(outnodes, outleaves, next, nextlevel);
		}
	}
}

//re-write levels in depth first order; leaves get their own file
void GSSTreeCreator::optimizeLevels()
{
	for(unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		nodes[i].switchToMap();
	}

	filesystem::path npath = dbpath / "nodes";
	ofstream outnodes(npath.string().c_str());

	filesystem::path lpath = dbpath / "leaves";
	ofstream outleaves(lpath.string().c_str());

	const GSSNodeCommon* root = (GSSNodeCommon*)nodes.back().map->get_address();

	optimizeLevelsR(outnodes, outleaves, root, nodes.size()-1);
}

//top down partition
void GSSLevelCreator::createNextLevel(DataViewer& data, ostream* nodefile, vector<file_index>& nodeindices,
		ostream* treefile, vector<file_index>& treeindices)
{
	if(data.size() == 0)
		return;

	TopDownPartitioner *thispart = partitioner->create(&data);

	packingSize = nodePack;
	if(data.isLeaf()) //making leaf nodes
		packingSize = leafPack;

	outNodes = nodefile;
	outTrees = treefile;
	nodeIndices = &nodeindices;
	treeIndices = &treeindices;
	//recursively partition
	createNextLevelR(thispart);
	delete thispart;
}

void GSSLevelCreator::createNextLevelR(TopDownPartitioner *P)
{
	if(P->size() == 0)
		return;

	if(P->size() <= packingSize)
	{
		//bottom up pack
		const DataViewer* data = P->getData();
		vector<Cluster> clusters;
		vector<unsigned> indices;
		P->extractIndicies(indices);
		packer->pack(data, indices, clusters);

		//each cluster creates a new node
		for(unsigned c = 0, nc = clusters.size(); c < nc; c++)
		{
			//write out node
			nodeIndices->push_back((file_index)outNodes->tellp());
			treeIndices->push_back((file_index)outTrees->tellp());
			if(data->isLeaf())
			{
				//the children are all single trees
				GSSLeaf::writeLeaf(data, clusters[c], *outNodes, *outTrees);
			}
			else
			{
				GSSInternalNode::writeNode(data, clusters[c], *outNodes, *outTrees);
			}
		}
	}
	else
	{
		vector<TopDownPartitioner*> parts;
		P->partition(parts);

		for(unsigned i = 0, n = parts.size(); i < n; i++)
		{
			createNextLevelR(parts[i]);
			delete parts[i];
		}
	}
}

