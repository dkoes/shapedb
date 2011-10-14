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
public:
	DataViewer();
	//these are file ind
	virtual const MappableOctTree* getMSV(file_index i) const = 0;
	virtual const MappableOctTree* getMIV(file_index i) const = 0;
};

class GSSTreeCreator
{
	WorkFile objects;
	WorkFile leaves;

	vector<WorkFile> nodes;

	filesystem::path dbpath;

	static void createNextLevel(DataViewer& data, const vector<file_index>& indices,
			ostream& out, vector<file_index>& nextindices);

public:
	GSSTreeCreator() {}
	~GSSTreeCreator() {}

	bool create(filesystem::path dir, Object::iterator itr, float dim, float res);
};

#endif /* GSSTREECREATOR_H_ */
