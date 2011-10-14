/*
 * GSSTreeCreator.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "GSSTreeCreator.h"

DataViewer::DataViewer()
{

}

class LeafViewer: public DataViewer
{
	const char *ptr;
public:
	LeafViewer(void *data): ptr((const char*)data)
	{

	}

	virtual const MappableOctTree* getMSV(file_index i) const
	{
		return NULL;
	}

	virtual const MappableOctTree* getMIV(file_index i) const
	{
		return NULL;
	}
};

class NodeViewer: public DataViewer
{
	const char *ptr;
public:
	NodeViewer(void *data): ptr((const char*)data)
	{

	}

	virtual const MappableOctTree* getMSV(file_index i) const
	{
		return NULL;
	}

	virtual const MappableOctTree* getMIV(file_index i) const
	{
		return NULL;
	}
};

WorkFile::WorkFile(const char *name): map(NULL)
{
	file = new ofstream(name);
	mapping = new file_mapping(name, read_only);
}

WorkFile::~WorkFile()
{
	if(file) delete file;
	if(mapping) delete mapping;
	if(map) delete map;
}

void WorkFile::set(const char *name)
{
	clear();
	file = new ofstream(name);
	mapping = new file_mapping(name, read_only);
}

void WorkFile::switchToMap()
{
	file->close();
	map = new mapped_region(*mapping, read_only);
}

//deallocate and reset
void WorkFile::clear()
{
	if(file) delete file;
	if(mapping) delete mapping;
	if(map) delete map;
	file = NULL;
	mapping = NULL;
	map = NULL;
}



//return true if successfull
bool GSSTreeCreator::create(filesystem::path dir, Object::iterator itr, float dim, float res)
{
	//create directory
	if (!filesystem::create_directory(dir))
	{
		cerr << "Unable to create database directory ";
		return false;
	}
	dbpath = dir;

	filesystem::path objfile = dbpath / "objs";
	filesystem::path leavesfile = dbpath / "leaves";

	//write out objects and trees
	objects.set(objfile.file_string().c_str());
	leaves.set(leavesfile.file_string().c_str());
	vector<file_index> indices;
	unsigned cnt = 0;
	for( ; itr ; ++itr)
	{
		const Object& obj = *itr;
		file_index pos = objects.file->tellp();
		obj.write(*objects.file);

		//leaf object
		indices.push_back((file_index)leaves.file->tellp());

		MappableOctTree *tree = MappableOctTree::create(dim, res, obj);
		GSSLeaf::writeLeaf(*leaves.file, pos, tree);
		delete tree;
		cnt++;
	}

	//partition leaves into bottom level
	//setup level
	stringstream levelname;
	levelname << "level" << nodes.size();
	filesystem::path nodesfile = dbpath / levelname.str();
	nodes.push_back(WorkFile(nodesfile.file_string().c_str()));

	//map the data
	leaves.switchToMap();
	LeafViewer leafdata(leaves.map->get_address());

	vector<file_index> nextindices; nextindices.reserve(indices.size()/2);
	createNextLevel(leafdata, indices, *nodes.back().file, nextindices);

	while(nextindices.size() > 1)
	{
		//map current level
		nodes.back().switchToMap();
		void *currentLevel = nodes.back().map->get_address();
		NodeViewer nodedata(currentLevel);

		//setup next level
		levelname.str("level");
		levelname << nodes.size();
		nodesfile = dbpath / levelname.str();
		nodes.push_back(WorkFile(nodesfile.file_string().c_str()));

		//shuffle indices
		swap(indices, nextindices);
		nextindices.clear();

		createNextLevel(nodedata, indices, *nodes.back().file, nextindices);
	}

	//output general info
	filesystem::path infoname = dbpath / "info";
	ofstream info(infoname.file_string().c_str());
	info << dim <<" " << res <<" " << nodes.size() << " " << cnt << "\n";

	return true;
}
