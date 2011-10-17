/*
 * GSSTreeCreator.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "GSSTreeCreator.h"

class TreeViewer: public DataViewer
{
public:
	TreeViewer(void *data): DataViewer(data)
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

	virtual bool isTree() const
	{
		return true;
	}
};

class NodeViewer: public DataViewer
{
	const char *ptr;
public:
	NodeViewer(void *data): DataViewer(data)
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

	virtual bool isTree() const
	{
		return false;
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
	trees.set(leavesfile.file_string().c_str());
	vector<file_index> indices;
	unsigned cnt = 0;
	for( ; itr ; ++itr)
	{
		const Object& obj = *itr;
		file_index pos = objects.file->tellp();
		obj.write(*objects.file);

		//leaf object
		indices.push_back((file_index)trees.file->tellp());

		MappableOctTree *tree = MappableOctTree::create(dim, res, obj);
		GSSTree::writeTree(*trees.file, pos, tree);
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
	trees.switchToMap();
	TreeViewer leafdata(trees.map->get_address());

	vector<file_index> nextindices; nextindices.reserve(indices.size()/2);
	leveler->createNextLevel(leafdata, indices, *nodes.back().file, nextindices);

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

		leveler->createNextLevel(nodedata, indices, *nodes.back().file, nextindices);
	}

	//output general info
	filesystem::path infoname = dbpath / "info";
	ofstream info(infoname.file_string().c_str());
	info << dim <<" " << res <<" " << nodes.size() << " " << cnt << "\n";

	return true;
}


//top down partition
void GSSLevelCreator::createNextLevel(DataViewer& data, const vector<file_index>& indices,
		ostream& o, vector<file_index>& next)
{
	if(indices.size() == 0)
		return;

	TopDownPartitioner *thispart = partitioner->create(&data, indices);

	packingSize = nodePack;
	if(data.isTree()) //making leaf nodes
		packingSize = leafPack;

	out = &o;
	nextindices = &next;
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
		vector<Packer::Cluster> clusters;
		vector<file_index> indices;
		P->extractIndicies(indices);
		packer->pack(P->getData(), indices, clusters);
		vector<const void*> children;
		for(unsigned c = 0, nc = clusters.size(); c < nc; c++)
		{
			children.clear();
			children.reserve(clusters[c].indices.size());
			for(unsigned i = 0, n = clusters[c].indices.size(); i < n; i++)
			{
				file_index index = clusters[c].indices[i];
				children.push_back(P->getData()->getAddress(index));
			}

		}

		nextindices->push_back((file_index)out->tellp());
		if(P->getData()->isTree())
		{
			//creating leaves from trees
			GSSLeaf::writeLeaf(*out, children);
		}
		else
		{
			GSSInternalNode::writeNode(*out, children);
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

