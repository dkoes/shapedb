/*
 * GSSTreeStructures.h
 *
 *  Created on: Oct 14, 2011
 *      Author: dkoes
 *
 *  This contains data types/structures that are shared between the creator and the searcher
 *
 *  These can be mapped into a file.
 */

#ifndef GSSTREESTRUCTURES_H_
#define GSSTREESTRUCTURES_H_

#include "MappableOctTree.h"
#include "GSSTypes.h"
#include <cassert>


struct GSSDoubleTree
{
	unsigned msvOffset;
	unsigned char data[]; //miv than msv

	//write out tree data in a manner consistent with structure layout
	static void writeTreeS(ostream& out, MappableOctTree *miv, MappableOctTree *msv)
	{
		unsigned p = miv->bytes();
		out.write((char*)&p, sizeof(unsigned));
		miv->write(out);
		msv->write(out);
	}

	const MappableOctTree* getMIV() const
	{
		return (const MappableOctTree*)data;
	}

	const MappableOctTree* getMSV() const
	{
		return (const MappableOctTree*)&data[msvOffset];
	}
};

//header of leaf and internal nodes
struct GSSNodeCommon
{
	bool isLeaf: 1;
	unsigned N: 31;
};

//a GSSLeaf only needs to store a single tree for each object, and the positions
//are within the object file, not the nodes file
class GSSLeaf
{
	struct Child
	{
		file_index object_pos;
		MappableOctTree tree;
	};

	GSSNodeCommon info;
	unsigned child_positions[]; //variable length, offset of trees after tree_positions
	unsigned char data[]; //convenient data ptr

public:
	static void writeLeaf(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees);

};

//a GSSInternalNode stores both the MIV and MSV for each subnode, and points to their locations
class GSSInternalNode
{
	struct Child
	{
		file_index node_pos;
		unsigned MSVindex; //where in data MSV starts (MIV at 0)
		unsigned char data[];
	};

	GSSNodeCommon info;
	unsigned child_positions[]; //offset of children after node_positions
	unsigned char data[]; //convienent data pointer

public:
	static void writeNode(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees);

};






#endif /* GSSTREESTRUCTURES_H_ */
