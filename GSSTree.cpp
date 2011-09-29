/*
 * GSSTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "GSSTree.h"
#include "Timer.h"
#include <math.h>
#include <iostream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/multi_array.hpp>
#include "CommandLine2/CommandLine.h"

using namespace std;

cl::opt<bool> ScanCheck("scancheck",cl::desc("Perform a full scan to check results"),cl::Hidden);
cl::opt<unsigned> KBestInsertion("kinsert",cl::desc("Evaluate k insertion points"),cl::value_desc("k"),cl::init(0));
cl::opt<bool> KBestReshuffle("reshuffle", cl::desc("Reshuffle using k insertion points on split"), cl::init(false));
cl::opt<unsigned> ReshuffleLimit("reshuffle-limit", cl::desc("Limit number of reshufflings"), cl::init(0));
cl::opt<bool> InsertionLoad("insert-load",cl::desc("Load be simply inserting"), cl::init(false));

cl::opt<unsigned> LeafPack("leaf-pack", cl::desc("Cutoff to trigger leaf packing"), cl::init(256));
cl::opt<unsigned> NodePack("node-pack", cl::desc("Cutoff to trigger leaf packing"), cl::init(0));

cl::opt<unsigned> LeafMerge("leaf-merge", cl::desc("Cutoff to merge small leaf clusters"), cl::init(4));
cl::opt<unsigned> NodeMerge("node-merge", cl::desc("Cutoff to merge small node clusters"), cl::init(4));

cl::opt<bool> HDist("hdist", cl::desc("Hausdorff distance as a splitting metric"), cl::init(false));

extern cl::opt<bool> Verbose;

const unsigned GSSTree::MaxSplit = 8; //number of children in each node

GSSTree::GSSNode::~GSSNode() { delete MSV; delete MIV; }


void GSSTree::setBoundingBox(array<float,6>& box)
{
	float dims[3] = {0,0,0};
	for(unsigned i = 0; i < 3; i++)
	{
		min[i] = box[2*i];
		dims[i] = box[2*i+1]-min[i];
	}

	//to make things easier, cubify using the longest dimension
	dim = max(dims[2],max(dims[0],dims[1]));
	//round all dimensions up to a power of 2; use fp
	dim = pow(2, ceil(log(dim)/log(2)));
	if(dim < 1)
		dim = 1;

	float div = dim;
	maxlevel = 0;
	while(div > maxres)
	{
		maxlevel++;
		div /=2;
	}
}

// perform any transformations necessary to conver the input coordinate system
//into the gss tree coordinate system
void GSSTree::transformMol(const vector<MolSphere>& m, vector<MolSphere>& ret)
{
	ret.clear();
	ret.reserve(m.size());

	for(unsigned i = 0, n = m.size(); i < n; i++)
	{
		ret.push_back(MolSphere(m[i].x-min[0], m[i].y-min[1], m[i].z-min[2],m[i].r));
	}
}

//add a single mol
void GSSTree::add(const vector<MolSphere>& m)
{
	//reposition mol
	vector<MolSphere> mol;
	transformMol(m, mol);
	OctTree *oct = octGen.newOctTree(dim, maxres, mol);

	//insert into the tree - first find the leaf node to insert into, then
	//perform the insertion

	float dist = HUGE_VAL;
	GSSLeafNode *leaf = NULL;

	vector<LeafDistPair> kbest;
	if(KBestInsertion > 0)
	{
		root->findInsertionPoints(oct, kbest, KBestInsertion);

		assert(kbest.size() > 0);
		//find the one that creates the smallest change in the MIV/MSV of the leaf
		unsigned besti = 0;
		float smallestChange = HUGE_VAL;
		for(unsigned i = 0, n = kbest.size(); i < n; i++)
		{
			float change = kbest[i].leaf->combinedVolumeChange(oct, oct);
			if(change < smallestChange)
			{
				smallestChange = change;
				besti = i;
			}
		}
		assert(isfinite(smallestChange));
		leaf = kbest[besti].leaf;

		if(!KBestReshuffle)
			kbest.clear();
	}
	else //single best
	{
		root->findInsertionPoint(oct, dist, leaf);
	}

	assert(leaf != NULL);
	leaf->insert(*this, oct, LeafData(m), kbest);
}

//initialize the index
GSSTree::PartitionData::PartitionData()
{

}

//compute shannon entropy of bit pattern distribution
// entropy = - Sum(p_i * log(p_i))
double GSSTree::Partitioner::shannonEntropy(unsigned *patternCnts, unsigned total)
{
	double res = 0;
	double tot = total;
	for(unsigned i = 0; i < 256; i++)
	{
		double pi = patternCnts[i]/tot;
		if(pi > 0)
		{
			res += pi*log(pi);
		}
	}

	return -res;
}

//enumerate all octants at this level
//TODO: be smarter - some can be ignored
void GSSTree::Partitioner::getOctants(unsigned level, bool splitMSV, vector< vector<unsigned> >& octants)
{
	//the number of octants grows exponentially with the level, for now
	//do the brute force approach, but eventually replace with a fraction ptr-based octtree
	//approach that traverse only the variable-presence octants
	octants.clear();
	octants.push_back( vector<unsigned>());
	for(unsigned l = 0; l < level; l++)
	{
		vector< vector<unsigned> > newoctants;
		newoctants.reserve(octants.size() * 8);
		for(unsigned i = 0; i < 8; i++)
		{
			for(unsigned j = 0, n = octants.size(); j < n; j++)
			{
				newoctants.push_back(octants[j]);
				newoctants.back().push_back(i);
			}
		}
		octants.swap(newoctants);
	}
}

// find an octant at the specified level that is "good"
// if splitMSV is true, look only at MSV patterns, otherwise MIV
// put result in octantcoord
bool GSSTree::Partitioner::findOctant(unsigned level, bool splitMSV, vector<unsigned>& octantcoord)
{
	//some ideas for "best" octant
	//-largest variance in volume
	//-most populated bit patterns
	//-shannon index of bit patterns  <----
	//-other diversity indices of bit patterns

	vector< vector<unsigned> > octants;
	getOctants(level, splitMSV, octants);

	double bestval = 0;
	unsigned besto = 0;
	for(unsigned o = 0, no = octants.size(); o < no; o++)
	{
		unsigned patternCnts[256] ={ 0, };
		unsigned total = 0;
		for (unsigned i = 0, n = tindex.size(); i < n; i++)
		{
			unsigned index = tindex[i];
			unsigned which = 0;
			if (splitMSV)
				which = partdata->getMSV(index)->getOctantPattern(octants[o], true);
			else
				which = partdata->getMIV(index)->getOctantPattern(octants[o], false);
			patternCnts[which]++;
			total++;
		}

		//various metrics are possible, probably just evaluate
		//shannon entropy and volume (faster?)
		double entropy = shannonEntropy(patternCnts, total);
		if(entropy > bestval)
		{
			bestval = entropy;
			besto = o;
		}
	}

	if(bestval == 0)
		return false; //all the same at this level

	octantcoord = octants[besto];
	return true;
}

//create sub partitions based on the octant bit pattern at the specified coord
void GSSTree::Partitioner::partitionOnOctant(const vector<unsigned>& octantcoord, bool splitMSV, vector<Partitioner>& parts)
{
	int pos[256];
	memset(pos, -1, sizeof(pos));

	parts.clear();
	parts.reserve(256);
	for(unsigned i = 0, n = tindex.size(); i < n; i++)
	{
		unsigned index = tindex[i];
		unsigned which = 0;
		if (splitMSV)
			which = partdata->getMSV(index)->getOctantPattern(octantcoord, true);
		else
			which = partdata->getMIV(index)->getOctantPattern(octantcoord, false);

		if(pos[which] < 0) //create new partition
		{
			pos[which] = parts.size();
			parts.push_back(Partitioner(partdata));
		}

		parts[pos[which]].addSingle(*this, i);
	}
}

//merge from into this, must have the same data
//does no checking for duplicates
void GSSTree::Partitioner::add(const Partitioner& from)
{
	assert(partdata == from.partdata);
	for(unsigned i = 0, n = from.size(); i < n; i++)
		tindex.push_back(from.tindex[i]);
}

//add a single bit of partition data from from, updating structures as necessary
void GSSTree::Partitioner::addSingle(const Partitioner& from, unsigned fromindex)
{
	tindex.push_back(from.tindex[fromindex]);
}

//pack the partition into clusters such that no cluster has more than max members
//presumably clusters should have no fewer than max/2 and there should be an
//effort made to keep the clusters dense, yet informative
//TODO: make efficient, evaluate clustering schemes
void GSSTree::Partitioner::packClusters(unsigned max, vector<Partitioner>& clusters)
{
	//compute pairwise distances between all members
	unsigned N = tindex.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);
	vector<vector<unsigned> > clusts; clusts.reserve(1+size() / max);
	for(unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const OctTree *imiv = partdata->getMIV(tindex[i]);
		const OctTree *imsv = partdata->getMSV(tindex[i]);
		for(unsigned j = 0; j < i; j++)
		{
			const OctTree *jmiv = partdata->getMIV(tindex[j]);
			const OctTree *jmsv = partdata->getMSV(tindex[j]);
			if(imiv == imsv && jmiv == jmsv) //leaf case
				distances[i][j] = distances[j][i] = leafDist(imiv, jmiv);
			else
				distances[i][j] = distances[j][i] = splitDist(imiv, imsv, jmiv,jmsv);
		}
		clusts.push_back(vector<unsigned>());
		clusts.back().push_back(i);
	}

	//iteratively merge all clusters until we run into size limit
	while(clusts.size() > 1)
	{
		//find two clusters with smallest complete linkage distance
		//(minimize the maximum distance)
		float mindist = HUGE_VAL;
		unsigned besti = 0;
		unsigned bestj = 0;
		for(unsigned i = 0, n = clusts.size(); i < n; i++)
		{
			for(unsigned j = 0; j < i; j++)
			{
				float maxdist = 0;
				//find the max distance between i and j
				for(unsigned I = 0, ni = clusts[i].size(); I < ni; I++)
				{
					for(unsigned J = 0, nj = clusts[j].size(); J < nj; J++)
					{
						float d = distances[clusts[i][I]][clusts[j][J]];
						if(d > maxdist)
						{
							maxdist = d;
						}
					}
				}
				if(maxdist < mindist && clusts[i].size() + clusts[j].size() <= max)
				{
					mindist = maxdist;
					besti = i;
					bestj = j;
				}
			}
		}

		if(besti == bestj)
			break; //couldn't merge any

		//move besti and bestj to end of vector
		unsigned last = clusts.size() -1;
		if(bestj == last)
		{
			swap(clusts[besti], clusts[last-1]);
		}
		else
		{
			swap(clusts[besti],clusts[last]);
			swap(clusts[bestj],clusts[last-1]);
		}

		//insert into second to last
		clusts[last-1].insert(clusts[last-1].end(), clusts[last].begin(), clusts[last].end());
		clusts.pop_back(); //remove
	}

	if(Verbose) cout << "  pack " << tindex.size() << " " << clusts.size() << "\n";

	//now create partitions
	clusters.clear();
	clusters.reserve(clusts.size());
	for(unsigned i = 0, n = clusts.size(); i < n; i++)
	{
		clusters.push_back(Partitioner(partdata));
		for(unsigned j = 0, m = clusts[i].size(); j < m; j++)
		{
			clusters.back().addSingle(*this, clusts[i][j]);
		}
	}
}

//create a leaf node from the contents of partitioner, does not check size
GSSTree::GSSLeafNode* GSSTree::leafFromPartition(Partitioner& partitioner, LeafPartitionData& leafdata)
{
	GSSLeafNode *newnode = new GSSLeafNode(octGen, dim, maxres);

	for(unsigned i = 0, n = partitioner.size(); i < n; i++)
	{
		newnode->trees.push_back(leafdata.getTree(partitioner.getDataIndex(i)));
		newnode->data.push_back(leafdata.getData(partitioner.getDataIndex(i)));
	}
	newnode->selfUpdate(); //update MIV/MSV

	return newnode;
}

//create an internal node from the contents of partitioner, does not check size
GSSTree::GSSInternalNode* GSSTree::nodeFromPartition(Partitioner& partitioner, NodePartitionData& nodedata)
{
	GSSInternalNode *newnode = new GSSInternalNode(octGen, dim, maxres);

	for(unsigned i = 0, n = partitioner.size(); i < n; i++)
	{
		newnode->children.push_back(nodedata.getNode(partitioner.getDataIndex(i)));
	}
	newnode->selfUpdate(); //update MIV/MSV

	return newnode;
}

/*
 * Take all the trees indexed by tindex, segregate them
 */
void GSSTree::partitionLeaves(Partitioner& partitioner, LeafPartitionData& leafdata, unsigned level, bool splitMSV, vector<GSSNode*>& nodes)
{
	if(partitioner.size() == 0)
		return;

	if(partitioner.size() <= MaxSplit)
	{
		//create a node
		nodes.push_back(leafFromPartition(partitioner, leafdata));
		return;
	}

	if(Verbose) cout << "leaves " << partitioner.size() << " l" << level << " s" << splitMSV << "\n";
	//identify octant to split on
	vector<unsigned> octantcoord;
	while(level <= maxlevel)
	{
		if(partitioner.findOctant(level, splitMSV, octantcoord))
			break;
		splitMSV = !splitMSV;
		if(partitioner.findOctant(level, splitMSV, octantcoord))
			break;
		level++;
	}

	//check for cases where we should just agglomerate
	if(level > maxlevel || partitioner.size() < LeafPack)
	{
		vector< Partitioner > clusters;
		partitioner.packClusters(MaxSplit, clusters);
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(leafFromPartition(clusters[i], leafdata));
		}
	}
	else
	{
		//first partition based on the occupancy pattern of the specified octant
		vector<Partitioner> parts;
		partitioner.partitionOnOctant(octantcoord, splitMSV, parts);

		//select which bins to recursively split and which to combine
		//into a single group to be clustered into leaves
		vector<unsigned> partitionMore;
		Partitioner mergedPartition(&leafdata);

		for(unsigned i = 0, n = parts.size(); i < n; i++)
		{
			if(parts[i].size() == 0)
				continue;
			if(parts[i].size() < LeafMerge)
			{
				mergedPartition.add(parts[i]);
			}
			else
			{
				partitionMore.push_back(i);
			}
		}

		//create leaves from small groups
		vector< Partitioner > clusters;
		mergedPartition.packClusters(MaxSplit, clusters);
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(leafFromPartition(clusters[i], leafdata));
		}

		//split on octants with large partitions
		for(unsigned i = 0, n = partitionMore.size(); i < n; i++)
		{
			partitionLeaves(parts[partitionMore[i]],leafdata, level, splitMSV, nodes);
		}
	}
}

//partition nodes to build another level
//very similar to partition leaves, but different data structures
//also, opportunaty for different algorithmic choices
void GSSTree::partitionNodes(Partitioner& partitioner, NodePartitionData& nodedata, unsigned level, bool splitMSV, vector<GSSNode*>& nodes)
{
	if(partitioner.size() == 0)
		return;
	if(partitioner.size() <= MaxSplit)
	{
		//create a node
		nodes.push_back(nodeFromPartition(partitioner, nodedata));
		return;
	}

	if(Verbose) cout << "nodes " << partitioner.size() << " l" << level << " s" << splitMSV << "\n";

	//identify octant to split on
	vector<unsigned> octantcoord;
	while(level <= maxlevel)
	{
		if(partitioner.findOctant(level, splitMSV, octantcoord))
			break;
		splitMSV = !splitMSV;
		if(partitioner.findOctant(level, splitMSV, octantcoord))
			break;
		level++;
	}

	//check for cases where we should just agglomerate
	if(level > maxlevel || partitioner.size() < NodePack)
	{
		vector< Partitioner > clusters;
		partitioner.packClusters(MaxSplit, clusters);
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(nodeFromPartition(clusters[i], nodedata));
		}
	}
	else
	{
		//first partition based on the occupancy pattern of the specified octant
		vector<Partitioner> parts;
		partitioner.partitionOnOctant(octantcoord, splitMSV, parts);

		//select which bins to recursively split and which to combine
		//into a single group to be clustered into leaves
		vector<unsigned> partitionMore;
		Partitioner mergedPartition(&nodedata);

		for(unsigned i = 0, n = parts.size(); i < n; i++)
		{
			if(parts[i].size() == 0)
				continue;
			if(parts[i].size() < NodeMerge)
			{
				mergedPartition.add(parts[i]);
			}
			else
			{
				partitionMore.push_back(i);
			}
		}

		//create leaves from small groups
		vector< Partitioner > clusters;
		mergedPartition.packClusters(MaxSplit, clusters);
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(nodeFromPartition(clusters[i], nodedata));
		}

		//split on octants with large partitions
		for(unsigned i = 0, n = partitionMore.size(); i < n; i++)
		{
			partitionNodes(parts[partitionMore[i]],nodedata, level, splitMSV, nodes);
		}
	}
}


void GSSTree::load(const vector<vector<MolSphere> >& mols)
{
	if (InsertionLoad)
	{
		//simply add molecules one at a time
		for (unsigned i = 0, n = mols.size(); i < n; i++)
		{
			add(mols[i]);
		}
	}
	else
	{
		//bulk loading

		//recursively partition trees by splitting on low resolution representations

		vector<OctTree*> trees;
		vector<LeafData> data;
		trees.reserve(mols.size());
		data.reserve(mols.size());
		for(unsigned i = 0, n = mols.size(); i < n; i++)
		{
			vector<MolSphere> mol;
			transformMol(mols[i], mol);

			trees.push_back(octGen.newOctTree(dim, maxres, mol));
			data.push_back(LeafData(mol));
		}

		LeafPartitionData leafdata(&trees, &data);
		Partitioner leafpartition(&leafdata);
		leafpartition.initFromData();
		vector<GSSNode*> nodes;
		vector<unsigned> octant;
		partitionLeaves(leafpartition, leafdata, 0, true, nodes);
		trees.clear(); //now stored in nodes
		data.clear();

		while(nodes.size() > 1)
		{
			NodePartitionData nodedata(&nodes);
			Partitioner nodepartition(&nodedata);
			nodepartition.initFromData();
			vector<GSSNode*> nextlevel;
			partitionNodes(nodepartition, nodedata, 0, true, nextlevel);
			swap(nextlevel, nodes);
		}
		assert(nodes.size() == 1);
		root = nodes[0];
	}
}


static unsigned leavesChecked = 0;
static unsigned fitChecks = 0;
//nearest neighbor search, return closest set of molspheres
void GSSTree::nn_search(const vector<MolSphere>& m, vector<MolSphere>& res)
{
	vector<MolSphere> mol;
	transformMol(m, mol);
	OctTree* tree = octGen.newOctTree(dim, maxres, mol);
	LeafData data;
	float dist = HUGE_VAL;
	leavesChecked = 0;
	fitChecks = 0;
	Timer t;
	root->findNearest(tree, dist, data);
	double elapsed = t.elapsed();
	if(Verbose) cout << "LeavesChecked " << leavesChecked << " / " << root->numLeaves() << "\tFitChecks " << fitChecks << " " << elapsed << "\n";
	swap(res, data.spheres);

	if(ScanCheck)
	{
		//check
		t.restart();
		leavesChecked = 0;
		fitChecks = 0;
		dist = HUGE_VAL;
		root->scanNearest(tree, dist, data);
		elapsed = t.elapsed();
		if(Verbose) cout << "LeavesScanned " << leavesChecked << " / " << root->numLeaves()
				<< "\tFitChecks " << fitChecks << " " << elapsed << "\n";

		if (data.spheres != res)
		{
			cout << "Find and Scan differ!\n";
		}
	}

	delete tree;
}

void GSSTree::tree_range_search(const OctTree* smallTree, const OctTree* bigTree, vector<vector<MolSphere>  >& res)
{
	leavesChecked = 0;
	fitChecks = 0;
	Timer t;
	res.clear();
	vector<LeafData> data;
	root->findTweeners(smallTree, bigTree, data);
	double elapsed = t.elapsed();
	if(Verbose) cout << "LeavesChecked "  << leavesChecked << " / " << root->numLeaves() << "\tFitChecks "
			<< fitChecks << " "<< elapsed << "\t" << data.size() << "\n";

	for(unsigned i = 0, n = data.size(); i < n; i++)
	{
		res.push_back(data[i].spheres);
	}

	if(ScanCheck)
	{
		//check
		t.restart();
		leavesChecked = 0;
		fitChecks = 0;
		data.clear();
		root->scanTweeners(smallTree, bigTree, data);
		elapsed = t.elapsed();
		if (Verbose) cout << "LeavesScanned " << leavesChecked << " / " << root->numLeaves()
								<< "\tFitChecks " << fitChecks << " " << elapsed << "\t" << data.size() << "\n";

		if (data.size() != res.size()) //lame but easy check
		{
			cout << "Find and Scan differ!\n";
		}
	}
}

void GSSTree::dc_search(const vector<MolSphere>& little, const vector<MolSphere>& big, vector<vector<MolSphere> >& res)
{
	vector<MolSphere> littleMol, bigMol;
	transformMol(little, littleMol);
	transformMol(big, bigMol);

	OctTree *smallTree = octGen.newOctTree(dim,maxres, littleMol);
	OctTree *bigTree = octGen.newOctTree(dim,maxres, bigMol);

	tree_range_search(smallTree, bigTree, res);

	delete smallTree;
	delete bigTree;
}

void GSSTree::inex_search(const vector<MolSphere>& inc, const vector<MolSphere>& exc, vector<vector<MolSphere> >& res)
{
	vector<MolSphere> littleMol, bigMol;
	transformMol(inc, littleMol);
	transformMol(exc, bigMol);

	OctTree* smallTree = octGen.newOctTree(dim,maxres, littleMol);
	OctTree* bigTree = octGen.newOctTree(dim,maxres, bigMol);
	bigTree->invert();

	tree_range_search(smallTree, bigTree, res);

	delete smallTree;
	delete bigTree;
}



//output leaf data
void GSSTree::LeafData::write(ostream& out) const
{
	unsigned n = spheres.size();
	out.write((char*)&n, sizeof(n));
	for(unsigned i = 0; i < n; i++)
	{
		out.write((char*)&spheres[i], sizeof(MolSphere));
	}
}

//and input
void GSSTree::LeafData::read(istream&  in)
{
	unsigned n = 0;
	in.read((char*)&n, sizeof(n));
	spheres.resize(n);
	for(unsigned i = 0; i < n; i++)
	{
		in.read((char*)&spheres[i], sizeof(MolSphere));
	}
}

//base class output
void GSSTree::GSSNode::write(ostream& out) const
{
	streampos ret = out.tellp();

	//write out kind of node
	const GSSInternalNode* intnode = dynamic_cast<const GSSInternalNode*>(this);
	bool isInternal = intnode != NULL;

	out.write((char*)&isInternal, sizeof(isInternal));
	out.write((char*)&res, sizeof(res));
	out.write((char*)&which, sizeof(which));
	//write out MIV/MSV
	MIV->write(out);
	MSV->write(out);
}


GSSTree::GSSNode* GSSTree::GSSNode::readCreate(const OctTreeFactory& octGen, istream& in, GSSInternalNode *parPtr)
{
	bool isInternal = false;
	in.read((char*)&isInternal, sizeof(isInternal));

	if(isInternal)
	{
		GSSInternalNode *ret = new GSSInternalNode(octGen);
		ret->read(octGen, in, parPtr);
		return ret;
	}
	else
	{
		GSSLeafNode *ret = new GSSLeafNode(octGen);
		ret->read(octGen, in, parPtr);
		return ret;
	}

}

void GSSTree::GSSNode::read(istream& in, GSSInternalNode *parPtr)
{
	parent = parPtr;
	in.read((char*)&res, sizeof(res));
	in.read((char*)&which, sizeof(which));
	MIV->read(in);
	MSV->read(in);
}

void GSSTree::GSSInternalNode::write(ostream& out) const
{
	GSSNode::write(out);

	unsigned n = children.size();
	out.write((char*)&n, sizeof(n));

	for(unsigned i = 0; i < n; i++)
	{
		children[i]->write(out);
	}
}

void GSSTree::GSSInternalNode::read(const OctTreeFactory& octGen, istream& in, GSSInternalNode *parPtr)
{
	GSSNode::read(in, parPtr);

	unsigned n = 0;
	in.read((char*)&n, sizeof(n));
	children.resize(n, NULL);
	for(unsigned i = 0; i < n; i++)
	{
		children[i] = GSSNode::readCreate(octGen, in, this);
	}
}

//dump leaf node
void GSSTree::GSSLeafNode::write(ostream& out) const
{
	GSSNode::write(out);

	unsigned n = trees.size();
	out.write((char*)&n, sizeof(n));
	for(unsigned i = 0; i < n; i++)
	{
		trees[i]->write(out);
	}

	for(unsigned i = 0; i < n; i++)
	{
		data[i].write(out);
	}
}

void GSSTree::GSSLeafNode::read(const OctTreeFactory& octGen, istream& in, GSSInternalNode *parPtr)
{
	GSSNode::read(in, parPtr);

	unsigned n = 0;
	in.read((char*)&n, sizeof(n));
	trees.resize(n, NULL);
	data.resize(n);

	for(unsigned i = 0; i < n; i++)
	{
		trees[i] = octGen.newOctTree();
		trees[i]->read(in);
	}

	for(unsigned i = 0; i < n; i++)
	{
		data[i].read(in);
	}
}


//dump tree in binary to be read into memory later
void GSSTree::write(ostream& out)
{
	out.write((char*)&maxres, sizeof(maxres));
	out.write((char*)min, sizeof(min));
	out.write((char*)&dim, sizeof(dim));
	octGen.write(out);
	root->write(out);
}

//read completely into memory
void GSSTree::read(istream& in)
{
	in.read((char*)&maxres, sizeof(maxres));
	in.read((char*)min, sizeof(min));
	in.read((char*)&dim, sizeof(dim));
	octGen.read(in);
	if(root) delete root;
	root = GSSNode::readCreate(octGen, in, NULL);
}

//evaluate the goodness of moving to from->to
//this is _expensive_ as we compute the change in volume
//for MIV/MSV of both leaves
float GSSTree::deltaFit(const GSSLeafNode *to, GSSLeafNode *from, unsigned t)
{
	//never remove the last tree from a node
	if(from->trees.size() <= 1)
		return 0;

	//compute MIV/MSV of from node without t
	OctTree* nMIV = octGen.newOctTree(dim, maxres);
	OctTree* nMSV = octGen.newOctTree(dim, maxres);
	nMIV->fill();

	for(unsigned i = 0, n = from->trees.size(); i < n; i++)
	{
		if(i != t)
		{
			nMIV->intersect(from->trees[i]);
			nMSV->unionWith(from->trees[i]);
		}
	}

	float fromMIV = nMIV->volume() - from->MIV->volume();
	float fromMSV = from->MSV->volume() - nMSV->volume();

	float fromCombined = fromMIV+fromMSV;
	float toCombined = to->combinedVolumeChange(from->trees[t], from->trees[t]);

	delete nMIV;
	delete nMSV;

	//only favorable if we tighten from more than we expand to
	return fromCombined - toCombined;
}

//insert data into leaf node, if there isn't enough room, split
//call update on parents if need-be
//takes ownership of the tree memory
void GSSTree::GSSLeafNode::insert(GSSTree& gTree, OctTree* tree, const LeafData& m, vector<LeafDistPair>& kbest)
{
	if(data.size() < MaxSplit)
	{
		bool needUpdate = false;

		assert(data.size() == trees.size());
		//easy case, just add
		trees.push_back(tree);
		data.push_back(LeafData(m));

		//update MIV/MSV
		needUpdate |= MIV->intersect(tree);
		needUpdate |= MSV->unionWith(tree);

		if(needUpdate && parent != NULL)
		{
			parent->update(gTree, which, NULL);
		}
	}
	else //must split
	{
		//first add tree
		trees.push_back(tree);
		data.push_back(LeafData(m));

		//so that we can have one split function, have concept of MIV and MSV
		//even at leaf
		vector<unsigned> split1;
		vector<unsigned> split2;
		split(trees, trees, split1, split2);
		GSSLeafNode *newnode = new GSSLeafNode(gTree.octGen, gTree.dim, gTree.maxres);

		//put split2 in new node
		newnode->MIV->fill();
		newnode->data.resize(split2.size());
		for(unsigned i = 0, n = split2.size(); i < n; i++)
		{
			unsigned index = split2[i];
			newnode->trees.push_back(trees[index]);
			swap(newnode->data[i],data[index]);

			//update MIV/MSV
			newnode->MSV->unionWith(trees[index]);
			newnode->MIV->intersect(trees[index]);
		}
		//now mogirfy to split1
		vector<OctTree*> tmptrees; tmptrees.reserve(MaxSplit);
		vector<LeafData> tmpdata(split1.size());

		MIV->fill();
		MSV->clear();
		for(unsigned i = 0, n = split1.size(); i < n; i++)
		{
			unsigned index = split1[i];
			tmptrees.push_back(trees[index]);
			swap(tmpdata[i],data[index]);

			MIV->intersect(trees[index]);
			MSV->unionWith(trees[index]);
		}

		swap(tmpdata,data);
		swap(tmptrees,trees);

		if(parent != NULL) parent->update(gTree, which, newnode);
		else if(newnode != NULL) //need new root
		{
			gTree.createRoot(this, newnode);
		}

		//if the kbest leaves were passed in, refine all the contents
		//of these leaves be evaluating whether or not each individual
		//tree would be more specific to one of the split nodes and moving it
		//if so - this is not cheap
		int cnt = ReshuffleLimit;
		if(kbest.size() > 0 && newnode != NULL)
		{
			for(unsigned i = 0, n = kbest.size(); i < n; i++)
			{
				if(cnt <= 0)
					break;
				GSSLeafNode *altleaf = kbest[i].leaf;
				if(altleaf == this)
					continue;

				//don't resplit, just stop if full
				if(trees.size() >= MaxSplit && newnode->trees.size() >= MaxSplit)
					break;

				bool changed1 = false, changed2 = false;
				for(int i = 0; i < (int)altleaf->trees.size(); i++)
				{
					float deltaFit1 = gTree.deltaFit(this, altleaf, i);
					float deltaFit2 = gTree.deltaFit(newnode, altleaf, i);

					//must be positive to be better fit
					if(deltaFit1 > 0 || deltaFit2 > 0)
					{
						if(deltaFit1 > deltaFit2 && trees.size() < MaxSplit)
						{
							//move to this node
							moveTreeFrom(altleaf, i);
							changed1 = true;
							i--; //has been deleted
							cnt--;
						}
						else if(newnode->trees.size() < MaxSplit)
						{
							//move to new node
							newnode->moveTreeFrom(altleaf, i);
							changed2 = true;
							i--;
							cnt--;
						}
					}
				}

				if(changed1 && parent != NULL) //update up the tree
					parent->fullUpdate();
				if(changed2 && newnode->parent != NULL)
					newnode->parent->fullUpdate();

			}
		}
	}

}

//add a child without any checks, just get the pointers right
void GSSTree::GSSInternalNode::addChild(GSSNode *child)
{
	child->parent = this;
	child->which = children.size();
	children.push_back(child);
}

//complete update all the way up the tree
void GSSTree::GSSNode::fullUpdate()
{
	selfUpdate();
	if(parent != NULL)
		parent->fullUpdate();
}

//update our MIV/MSV from scratch (necessary when removing)
//do not recurse upwards
void GSSTree::GSSLeafNode::selfUpdate()
{
	MIV->fill();
	MSV->clear();

	for(unsigned i = 0, n = trees.size(); i < n; i++)
	{
		MIV->intersect(trees[i]);
		MSV->unionWith(trees[i]);
	}
}

//update our MIV/MSV from scratch (necessary when removing)
//do not recurse upwards
void GSSTree::GSSInternalNode::selfUpdate()
{
	MIV->fill();
	MSV->clear();

	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		MIV->intersect(children[i]->MIV);
		MSV->unionWith(children[i]->MSV);
	}
}

//recurse UP the tree updating MIV/MSV and splitting as necessary
void GSSTree::GSSInternalNode::update(GSSTree& tree, unsigned whichChild, GSSNode *newnode)
{
	//first update MIV/MSV
	assert(whichChild < children.size());

	bool needUpdate = false;
	needUpdate |= MIV->intersect(children[whichChild]->MIV);
	needUpdate |= MSV->unionWith(children[whichChild]->MSV);

	if(newnode != NULL)
	{
		//add node to this level
		if(children.size() < MaxSplit)
		{
			//easy case
			needUpdate |= MIV->intersect(newnode->MIV);
			needUpdate |= MSV->unionWith(newnode->MSV);
			addChild(newnode);
			newnode = NULL;
		}
		else //must split
		{
			//setup split structures
			vector<OctTree*> MIVs, MSVs;
			MIVs.reserve(MaxSplit);
			MSVs.reserve(MaxSplit);

			//add in new node briefly so it's part of the split
			children.push_back(newnode);

			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				MIVs.push_back(children[i]->MIV);
				MSVs.push_back(children[i]->MSV);
			}

			vector<unsigned> split1, split2;

			//do the split
			split(MIVs, MSVs, split1, split2);

			//save oldchildren
			vector<GSSNode*> oldchildren;
			oldchildren.reserve(MaxSplit);
			swap(oldchildren, children);

			//mogify this node to hold split1
			MIV->fill();
			MSV->clear();
			for(unsigned i = 0, n = split1.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split1[i]];
				addChild(child);
				MIV->intersect(child->MIV);
				MSV->unionWith(child->MSV);
			}

			//same thing for the new node and split2
			GSSInternalNode *newinode = new GSSInternalNode(tree.octGen,tree.dim,tree.maxres);
			newinode->MIV->fill();
			newinode->MSV->clear();
			for(unsigned i = 0, n = split2.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split2[i]];
				newinode->addChild(child);
				newinode->MIV->intersect(child->MIV);
				newinode->MSV->unionWith(child->MSV);
			}
			newnode = newinode;
		}
	}

	if(parent != NULL)
	{
		parent->update(tree, which, newnode);
	}
	else if(newnode != NULL && parent == NULL) //need new root
	{
		tree.createRoot(this, newnode);
	}

}

//create a new root using left and right, which are assumed to
//be splits of the current root
void GSSTree::createRoot(GSSNode *left, GSSNode *right)
{
	assert(left == root);
	GSSInternalNode* newroot = new GSSInternalNode(octGen,dim, maxres);
	//update trees
	*newroot->MIV = *left->MIV;
	newroot->MIV->intersect(right->MIV);

	*newroot->MSV = *left->MSV;
	newroot->MSV->unionWith(right->MSV);

	//add two children
	newroot->addChild(left);
	newroot->addChild(right);

	root = newroot;
}


//return true if the object(s) represented by MIV/MSV might fit in between min and max
bool GSSTree::fitsInbetween(const OctTree *MIV, const OctTree *MSV, const OctTree *min, const OctTree *max)
{
	fitChecks++;
	//TODO: more efficient
	//the MSV must completely enclose min
	//cout << MSV->volume() << " , " << min->volume() << "," << MSV->unionVolume(min) << "\n";
	if(MSV->unionVolume(min) != MSV->volume())
		return false;
	//MIV must be completely enclosed by max
	if(MIV->intersectVolume(max) != MIV->volume())
		return false;

	return true;
}

//return a distance to a single leaf, should be compatible with search Dist
float GSSTree::leafDist(const OctTree* obj, const OctTree *leaf)
{
	if(HDist)
	{
		return max(obj->hausdorffDistance(leaf), leaf->hausdorffDistance(obj));
	}
	else
	{
		return 1 - obj->intersectVolume(leaf)/obj->unionVolume(leaf);
	}
}

//return a "distance" between obj and MIV/MSV; the lower the distance
//the higher the more likely a node should be searched
//min and max should bookend the ultimate leaf distances
float GSSTree::searchDist(const OctTree* obj, const OctTree *MIV, const OctTree *MSV, float& min, float& max)
{
	min = 1 - obj->intersectVolume(MSV)/obj->unionVolume(MIV);
	max = 1 - obj->intersectVolume(MIV)/obj->unionVolume(MSV);

	return min+max;
}


/*
 * return a "distance" between approximations
 * for now, just the sum of the volume different of the MIVs and of the MSVs
 */
float GSSTree::splitDist(const OctTree* leftMIV, const OctTree* leftMSV, const OctTree* rightMIV, const OctTree* rightMSV)
{
	if(HDist)
	{
		float d1 = max(leftMIV->hausdorffDistance(rightMIV),rightMIV->hausdorffDistance(leftMIV));
		float d2 = max(leftMSV->hausdorffDistance(rightMSV),rightMSV->hausdorffDistance(leftMSV));

		return d1+d2;
	}
	else
	{
		float mind = 1 - leftMIV->intersectVolume(rightMIV)
				/ leftMIV->unionVolume(rightMIV);
		float maxd = 1 - leftMSV->intersectVolume(rightMSV)
				/ leftMSV->unionVolume(rightMSV);

		return mind + maxd;
	}
}


/* Decide how to split a node.  For now, use the canonical quadratic split
 * where you find the two most "distant" nodes, and then greedily partition
 * using those.
 */
void GSSTree::split(const vector<OctTree*>& MIV, const vector<OctTree*>& MSV,
		vector<unsigned>& s1, vector<unsigned>& s2)
{
	assert(MIV.size() == MSV.size());
	s1.clear();
	s2.clear();
	unsigned besti = 0;
	unsigned bestj = 0;
	float dist = 0;

	unsigned N = MIV.size();
	s1.reserve(N);
	s2.reserve(N);
	float distances[N][N];

	for(unsigned i = 0; i < N; i++)
	{
		for(unsigned j = i+1; j < N; j++)
		{
			float d = splitDist(MIV[i],MSV[i],MIV[j],MSV[j]);
			distances[i][j] = d;
			distances[j][i] = d;
			if(d > dist)
			{
				dist = d;
				besti = i;
				bestj = j;
			}
		}
		distances[i][i] = 0;
	}
	assert(besti != bestj);

	//seed with besti and bestj
	OctTree *iMIV = MIV[besti]->clone();
	OctTree *iMSV = MSV[besti]->clone();
	s1.push_back(besti);

	OctTree *jMIV = MIV[bestj]->clone();
	OctTree *jMSV = MSV[bestj]->clone();
	s2.push_back(bestj);

	for(unsigned i = 0; i < N; i++)
	{
		if(i != besti && i != bestj)
		{
			//choose based on distance to cumulative MIV/MSV
			//unless we've already filled one split
			if(s1.size() > N/2+1)
			{
				s2.push_back(i);
			}
			else if(s2.size() > N/2+1)
			{
				s1.push_back(i);
			}
			else
			{
				float di = splitDist(iMIV,iMSV,MIV[i],MSV[i]);
				float dj = splitDist(jMIV, jMSV, MIV[i],MSV[i]);

				if(di < dj)
				{
					s1.push_back(i);
					iMIV->intersect(MIV[i]);
					iMSV->unionWith(MSV[i]);
				}
				else if(dj < di)
				{
					s2.push_back(i);
					jMIV->intersect(MIV[i]);
					jMSV->unionWith(MSV[i]);
				}
				else //equal, seems unlikely
				{
					if(s1.size() < s2.size())
					{
						s1.push_back(i);
						iMIV->intersect(MIV[i]);
						iMSV->unionWith(MSV[i]);
					}
					else
					{
						s2.push_back(i);
						jMIV->intersect(MIV[i]);
						jMSV->unionWith(MSV[i]);
					}
				}

			}
		}
	}


	delete iMIV;
	delete iMSV;
	delete jMIV;
	delete jMSV;
}

//how much will the current MIV and MSV change after merging in miv and msv?
float GSSTree::GSSNode::combinedVolumeChange(const OctTree *miv, const OctTree *msv) const
{
	float deltamiv = MIV->volume() - MIV->intersectVolume(miv);
	float deltamsv = MSV->unionVolume(msv) - MSV->volume();
	return deltamiv + deltamsv;
}


//examine every object in the leaf and see if it has a better distance than distance,
//if so, store data
void GSSTree::GSSLeafNode::scanNearest(const OctTree* tree, float& distance, LeafData& d)
{
	findNearest(tree, distance, d); //this is just the same
}


//find the object with the best distance in this node, if it's better than
//the passed distance, update
void GSSTree::GSSLeafNode::findNearest(const OctTree* tree, float& distance, LeafData& d)
{
	leavesChecked++;
	for(unsigned i = 0, n = trees.size(); i < n; i++)
	{
		float dist = leafDist(tree, trees[i]);
		if(dist < distance)
		{
			distance = dist;
			d = data[i];
		}
	}
}

//identify all objects that are exactly in between min and max
void GSSTree::GSSLeafNode::scanTweeners(const OctTree* min, const OctTree* max, vector<LeafData>& res)
{
	findTweeners(min, max, res);
}

static unsigned foundLeafCnt = 0;
//identify all objects that are exactly in between min and max
void GSSTree::GSSLeafNode::findTweeners(const OctTree* min, const OctTree* max, vector<LeafData>& res)
{
	leavesChecked++;
	bool found = false;
	for(unsigned i = 0, n = trees.size(); i < n; i++)
	{
		if(fitsInbetween(trees[i], trees[i],min,max))
		{
			found = true;
			res.push_back(data[i]);
		}
	}

	if(0 && found)
	{
		//output stuff
		cout << "FOUNDLEAF " << foundLeafCnt << "\n";
		stringstream name;
		name << "LEAF" << foundLeafCnt << ".xyz";
		ofstream out(name.str().c_str());

		for(unsigned i = 0, n = trees.size(); i < n; i++)
		{
			bool good = fitsInbetween(trees[i], trees[i],min,max);
			cout << i << " " << good << " " << splitDist(trees[i], trees[i], min, max) << " |";
			for(unsigned j = 0; j < n; j++)
			{
				cout << " " << splitDist(trees[j],trees[j],trees[i],trees[i]);
			}
			cout << "\n";

			//dump mol
			out << data[i].spheres.size();
			out << "\nLEAF" << foundLeafCnt << " " << i << "\n";
			for (unsigned j = 0, m = data[i].spheres.size(); j < m; j++)
			{
				out << "C " << data[i].spheres[j].x << " " << data[i].spheres[j].y << " "
						<< data[i].spheres[j].z << "\n";
			}
		}
		cout << "----\n";
		foundLeafCnt++;

	}
}


//if this leaf is a more appropriate place for tree, return self
void GSSTree::GSSLeafNode::findInsertionPoint(const OctTree* tree, float& distance, GSSLeafNode*& leaf)
{
	if(trees.size() == 0)
	{
		distance = 0;
		leaf = this;
	}
	else
	{
		float min, max;
		float dist = searchDist(tree, MIV, MSV, min, max);
		if(dist < distance)
		{
			distance = dist;
			leaf = this;
		}
	}
}

//add this leaf if it's distance is better and truncate kbest if necessary
void GSSTree::GSSLeafNode::findInsertionPoints(const OctTree* tree, vector<LeafDistPair>& kbest, unsigned k)
{
	float dist = 0;
	if(trees.size() > 0)
	{
		float min, max;
		dist = searchDist(tree, MIV, MSV, min, max);
	}

	if(kbest.size() == 0 || kbest.back().distance > dist)
	{
		LeafDistPair pair(this, dist);
		kbest.insert(lower_bound(kbest.begin(), kbest.end(), pair), pair);
		if(kbest.size() > k) //truncate
			kbest.resize(k);
	}

}

//move the t'th tree of from into this, removing from from
//updates the local MIV/MSV, but does not recurse up the tree
void GSSTree::GSSLeafNode::moveTreeFrom(GSSLeafNode* from, unsigned t)
{
	assert(t < from->trees.size());

	swap(from->trees[t],from->trees.back());
	swap(from->data[t], from->data.back());

	trees.push_back(from->trees.back());
	data.push_back(from->data.back());

	from->trees.pop_back();
	from->data.pop_back();

	from->selfUpdate();

	//update MIV/MSV
	MIV->intersect(trees.back());
	MSV->unionWith(trees.back());
}



struct ScoreIndex
{
	unsigned index;
	float score;
	float min;

	ScoreIndex() {}
	ScoreIndex(unsigned i, float s, float m): index(i), score(s), min(m) {}

	bool operator<(const ScoreIndex& si) const
	{
		return score < si.score;
	}
};


//identify all objects that are exactly in between min and max
//brute force scan for debugging
void GSSTree::GSSInternalNode::scanTweeners(const OctTree* min, const OctTree* max, vector<LeafData>& res)
{
	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		children[i]->scanTweeners(min, max, res);
	}
}

//identify all objects that are exactly in between min and max
//filter out children that can't possibly match
void GSSTree::GSSInternalNode::findTweeners(const OctTree* min, const OctTree* max, vector<LeafData>& res)
{
	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		if(fitsInbetween(children[i]->MIV, children[i]->MSV, min, max))
		{
			children[i]->findTweeners(min, max, res);
		}
	}
}

//scan - ignore any bounding on the search
void GSSTree::GSSInternalNode::scanNearest(const OctTree* tree, float& distance, LeafData& data)
{
	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		children[i]->scanNearest(tree, distance, data);
	}
}

//explore children to find closest value
void GSSTree::GSSInternalNode::findNearest(const OctTree* tree, float& distance, LeafData& data)
{
	vector<ScoreIndex> SIs;
	SIs.reserve(children.size());

	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		float min = 0, max = 0;
		float score = searchDist(tree, children[i]->MIV, children[i]->MSV, min, max);
		if(min < distance) //there's hope
		{
			SIs.push_back(ScoreIndex(i,score,min));
		}
	}

	//explore in priority order
	sort(SIs.begin(), SIs.end());

	for(unsigned i = 0, n = SIs.size(); i < n; i++)
	{
		if(SIs[i].min < distance)
		{
			children[SIs[i].index]->findNearest(tree, distance, data);
		}
	}
}

//look for a good place for the data indexed by tree
void GSSTree::GSSInternalNode::findInsertionPoint(const OctTree* tree, float& distance, GSSLeafNode*& leaf)
{
	float min, max;
	float dist = searchDist(tree, MIV, MSV, min, max);
	if(dist < distance)
	{
		//if this tree might contain something better, look at all children -> really should sort this
		for(unsigned i = 0, n = children.size(); i < n; i++)
		{
			children[i]->findInsertionPoint(tree, distance, leaf);
		}
	}
}

void GSSTree::GSSInternalNode::findInsertionPoints(const OctTree* tree, vector<LeafDistPair>& kbest, unsigned k)
{
	float min, max;
	float dist = searchDist(tree, MIV, MSV, min, max);
	float bestdist = HUGE_VAL;
	if(kbest.size() == k && k > 0)
	{
		bestdist = kbest.back().distance;
	}
	if(dist < bestdist)
	{
		//if this tree might contain something better, look at all children -> really should sort this
		for(unsigned i = 0, n = children.size(); i < n; i++)
		{
			children[i]->findInsertionPoints(tree, kbest, k);
		}
	}
}

//collect aggregate statistics for density of nodes and leaves
//returns depth
unsigned GSSTree::GSSInternalNode::getStats(Stats& leaves, Stats& nodes) const
{
	unsigned sz = children.size();
	if(sz > nodes.max)
		nodes.max = sz;
	if(sz < nodes.min)
		nodes.min = sz;
	nodes.total += sz;
	nodes.cnt++;
	if(sz == 1) nodes.singletonCnt++;
	unsigned ret = 0;
	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		ret = max(ret,children[i]->getStats(leaves,nodes));
	}
	return ret+1;
}

unsigned GSSTree::GSSLeafNode::getStats(Stats& leaves, Stats& nodes) const
{
	unsigned sz = trees.size();
	if(sz > leaves.max)
		leaves.max = sz;
	if(sz < leaves.min)
		leaves.min = sz;
	leaves.total += sz;
	leaves.cnt++;
	if(sz == 1) leaves.singletonCnt++;
	return 0;
}

//print out aggregrate statistics for density of nodes/leaves
void GSSTree::printStats() const
{
	Stats leaves, nodes;
	unsigned depth = root->getStats(leaves,nodes);

	printf("%dDepth %fAve | Leaves (%d): %fAve %dMin %dMax %dSingle| Nodes (%d): %fAve %dMin %dMax %dSingle\n",
			depth,(leaves.total+nodes.total)/(double)(leaves.cnt+nodes.cnt),
			leaves.cnt, leaves.total/(double)leaves.cnt, leaves.min, leaves.max, leaves.singletonCnt,
			nodes.cnt, nodes.total/(double)nodes.cnt, nodes.min, nodes.max, nodes.singletonCnt);
}



