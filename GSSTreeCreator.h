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

#include <iostream>
#include <fstream>
#include <vector>

#include "GSSTypes.h"
#include "GSSTreeStructures.h"
#include "TopDownPartitioner.h"
#include "packers/Packer.h"

using namespace boost;
using namespace std;

#include "WorkFile.h"
#include "Molecule.h"

typedef Molecule Object; //eventually template this

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
	ostream *outNodes;
	ostream *outTrees;
	vector<file_index> *nodeIndices;
	vector<file_index> *treeIndices;
	virtual void createNextLevelR(TopDownPartitioner *P);
public:
	GSSLevelCreator(const TopDownPartitioner * part, const Packer *pack,
			unsigned np, unsigned lp) :
			partitioner(part), packer(pack), nodePack(np), leafPack(lp)
	{
	}

	virtual ~GSSLevelCreator()
	{
	}

	virtual void createNextLevel(DataViewer& data, ostream* nodefile,
			vector<file_index>& nodeindices, ostream* treefile,
			vector<file_index>& treeindices);

	unsigned getPack() const
	{
		return packer->getPack();
	}
};

class GSSTreeCreator
{
	WorkFile objects;

	vector<WorkFile> nodes;

	filesystem::path dbpath;

	GSSLevelCreator *leveler;

	float dimension;
	float resolution;
	unsigned superNodeDepth;
	//some bookkeeping for analysis purposes
	unsigned numNodes;
	unsigned numLeaves;
	vector<unsigned> nodeContentDistribution;
	vector<unsigned> leafContentDistribution;
	file_index optimizeLevelsR(ostream& outnodes, ostream& outleaves,
			const GSSNodeCommon *n, unsigned level, file_index& lstart,
			file_index& lend);
	void optimizeLevels();

	void getNodesForSuperNode(const GSSInternalNode* root,
			vector<GSSInternalNode*>& newroots, unsigned curlevel,
			unsigned stoplevel);

public:
	GSSTreeCreator(GSSLevelCreator *l, unsigned sdepth=3) :
			leveler(l), dimension(0), resolution(0), superNodeDepth(sdepth)
	{
	}
	~GSSTreeCreator()
	{
	}

	bool create(filesystem::path dir, Object::iterator& itr, float dim,
			float res);

	void printStats(ostream& out) const;
};

#endif /* GSSTREECREATOR_H_ */
