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

struct GSSLeaf
{
	file_index object_pos;
	MappableOctTree tree;

	//write out leaf data in a manner consistent with structure layout
	static void writeLeaf(ostream& out, file_index pos, MappableOctTree *tree)
	{
		out.write((char*)&pos, sizeof(file_index));
		tree->write(out);
	}
};



#endif /* GSSTREESTRUCTURES_H_ */
