/*
 * GSSTreeCreator.h
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 *
 *      This class creates a GSS tree on disk.  It assumes the input
 *      does not fit in memory and behaves accordingly.
 *
 *      It takes an iterator over the input data, which must support
 *      an intersection (with a cube) method and a write to file method.
 *      This data is then converted to oct tree representations.  The
 *      data gets written to an indexed file while the oct trees get written
 *      to another file (storing indices to the object data).
 *
 *      The tree file is then clustered to create leaves which are
 *      appended to a leaf file as they are created.
 *
 *      The leaf nodes are similarily clustered to create a level of
 *      nodes, which are written to their own file.  These nodes are
 *      clustered into the next level's file and so on.
 *
 *      Clustering involves a top-down O(n) partitioning that splits the
 *      data until the a set small enough for an O(n^2) bottom-up packing
 *      to be performed.  Clusters are always packed to contain at least 2
 *      entries (this may be relaxed for leaves).
 *
 *      Once the final level is created, the internal nodes are all laid out
 *      in a separate file in depth-first order.  The leaves are ordered
 *      sequentially in a separate file.
 */

#ifndef GSSTREECREATOR_H_
#define GSSTREECREATOR_H_

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <vector>

#include "GSSTreeStructures.h"

using namespace boost;
using namespace std;
using namespace boost::interprocess;

#include "Molecule.h"
typedef Molecule Object;

//store infor for files that are being created and will be memory mapped
struct WorkFile
{
	//none of these has a real copy constructor
	ofstream *file;
	mapped_region *map;
	file_mapping *mapping;

	WorkFile(): file(NULL), map(NULL), mapping(NULL) {}
	WorkFile(const char *name);
	~WorkFile();

	void switchToMap();
	void set(const char *name);
	void clear();
};

//a wrapper that can view single tree leaves the same as internal nodes
class DataViewer
{
	const char *ptr;
public:
	DataViewer(void *data): ptr((const char*)data) {}
	virtual ~DataViewer() {}
	//these are file ind
	virtual const MappableOctTree* getMSV(file_index i) const = 0;
	virtual const MappableOctTree* getMIV(file_index i) const = 0;
	virtual bool isTree() const  = 0;

	const void* getAddress(file_index i) const { return ptr+i; }
};


//abstract class for a top down partitioner, is expected to run in linear time
//can maintain state as the partitions are refined
class TopDownPartitioner
{
protected:
	const DataViewer *data;
	vector<file_index> indices;
public:
	TopDownPartitioner() {}
	virtual ~TopDownPartitioner() {}
	virtual TopDownPartitioner* create(const DataViewer* dv, const vector<file_index>& indices) const = 0;
	virtual void partition(vector<TopDownPartitioner*>& parts) = 0;

	virtual unsigned size() const { return indices.size(); }
	virtual const DataViewer * getData() const { return data; }
	virtual void extractIndicies(vector<file_index>& ind) { swap(ind, indices); }
};

//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
public:
	struct Cluster
	{
		vector<file_index> indices;
		MappableOctTree *MIV;
		MappableOctTree *MSV;
	};

	Packer() {}
	virtual ~Packer() {}
	virtual void pack(const DataViewer* dv, const vector<file_index>& indices, vector<Cluster>& clusters) const;
};

//class for creating levels, follows the CM-tree bulk loading algorithm,
//but can be overridden to implement any arbitrary algorithm
class GSSLevelCreator
{
protected:
	const TopDownPartitioner *partitioner;
	const Packer *packer;

	//configuration settings
	unsigned nodePack;
	unsigned leafPack;


	//class vars used by nextlevelR
	unsigned packingSize;
	ostream *out;
	vector<file_index> *nextindices;
	virtual void createNextLevelR(TopDownPartitioner *P);
public:
	GSSLevelCreator(const TopDownPartitioner * part, const Packer *pack, unsigned np, unsigned lp):
		partitioner(part), packer(pack), nodePack(np), leafPack(lp) {}

	virtual ~GSSLevelCreator() {}

	virtual void createNextLevel(DataViewer& data, const vector<file_index>& indices,
			ostream& out, vector<file_index>& nextindices);
};

class GSSTreeCreator
{
	WorkFile objects;
	WorkFile trees;

	vector<WorkFile> nodes;

	filesystem::path dbpath;

	GSSLevelCreator *leveler;

public:
	GSSTreeCreator(GSSLevelCreator *l):leveler(l) {}
	~GSSTreeCreator() {}

	bool create(filesystem::path dir, Object::iterator itr, float dim, float res);
};

#endif /* GSSTREECREATOR_H_ */
