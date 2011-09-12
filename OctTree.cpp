/*
 * OctTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "OctTree.h"

//perform deep recursive copy
void OctTree::copy(OctNode *dst, const OctNode *src)
{
	dst->type = src->type;
	for(unsigned i = 0; i < 8; i++)
	{
		if(src->children[i] != NULL)
		{
			dst->children[i] = new OctNode();
			copy(dst->children[i], src->children[i]);
		}
	}
}

//deep copy constructor
OctTree::OctTree(const OctTree& rhs)
{
	resolution = rhs.resolution;
	dimension = rhs.dimension;

	copy(&root, &rhs.root);
}


//recursive creation of oct tree
void OctTree::create(OctNode *node, const Cube& cube, const vector<MolSphere>& mol)
{
	//does the mol overlap with this cube?
	vector<MolSphere> pruned;
	pruned.reserve(mol.size());
	for(unsigned i = 0, n = mol.size(); i < n; i++)
	{
		if(cube.intersectsSphere(mol[i]))
		{
			pruned.push_back(mol[i]);
		}
	}

	if(pruned.size() == 0)
	{
		//no overlap, all done
		node->type = Empty;
	}
	else if(cube.getDimension() <= resolution) //consider it full
	{
		node->type = Full;
	}
	else //subdivide into children
	{
		node->type = Children;
		unsigned fullcnt = 0;
		for(unsigned i = 0; i < 8; i++)
		{
			node->children[i] = new OctNode();
			Cube newc = cube.getOctant(i);
			create(node->children[i], newc, pruned);
			if(node->children[i]->type == Full)
				fullcnt++;
		}

		//are all the children full? then truncate and mark node as full
		if(fullcnt == 8)
		{
			node->deleteChildren();
			node->type = Full;
		}
	}
}
