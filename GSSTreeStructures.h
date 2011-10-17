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

typedef unsigned long file_index;

struct GSSTree
{
	file_index object_pos;
	MappableOctTree tree;

	//write out tree data in a manner consistent with structure layout
	static void writeTree(ostream& out, file_index pos, MappableOctTree *tree)
	{
		out.write((char*)&pos, sizeof(file_index));
		tree->write(out);
	}
};

//a GSSLeaf only needs to store a single tree for each object, and the positions
//are within the object file, not the nodes file
struct GSSLeaf
{
	struct Child
	{
		file_index object_pos;
		MappableOctTree tree;
	};

	bool isLeaf: 1;
	file_index object_pos: 63;
	unsigned N; //number of objects
	unsigned tree_positions[]; //variable length, offset of trees after tree_positions

	static void writeLeaf(ostream& out, const vector<const void*>& trees);
};

//a GSSInternalNode stores both the MIV and MSV for each subnode, and points to their locations
struct GSSInternalNode
{
	struct Child
	{
		file_index node_pos;
		unsigned MSVindex; //where in data MSV starts (MIV at 0)
		unsigned char data[];
	};

	bool isLeaf: 1;
	file_index object_pos: 63;
	unsigned N; //number of objects
	unsigned node_positions[]; //offset of children after node_positions

	static void writeLeaf(ostream& out, const vector<const void*>& nodes);

};

#endif /* GSSTREESTRUCTURES_H_ */
