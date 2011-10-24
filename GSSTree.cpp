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
#include <queue>
#include "CommandLine2/CommandLine.h"


using namespace std;

extern cl::opt<bool> ScanCheck;
cl::opt<unsigned> KBestInsertion("kinsert",
		cl::desc("Evaluate k insertion points"), cl::value_desc("k"),
		cl::init(0));
cl::opt<bool> KBestReshuffle("reshuffle",
		cl::desc("Reshuffle using k insertion points on split"),
		cl::init(false));
cl::opt<unsigned> ReshuffleLimit("reshuffle-limit",
		cl::desc("Limit number of reshufflings"), cl::init(0));
cl::opt<bool> InsertionLoad("insert-load",
		cl::desc("Load be simply inserting"), cl::init(false));

extern cl::opt<unsigned> LeafPack;
extern cl::opt<unsigned> NodePack;

cl::opt<unsigned> LeafMerge("leaf-merge",
		cl::desc("Cutoff to merge small leaf clusters"), cl::init(4));
cl::opt<unsigned> NodeMerge("node-merge",
		cl::desc("Cutoff to merge small node clusters"), cl::init(4));

cl::opt<bool> PackPartitions("pack-partitions", cl::desc("Group partitions for more even split"), cl::init(false));
cl::opt<bool> AdaptivePacking("adaptive-pack", cl::desc("Dynamically set packing amount"), cl::init(false));

enum SplitMetricEnum
{
	AverageRelVolume, AverageAbsVolume, AverageHausdorff, SeparationVolume,SharedOverlap
};

cl::opt<SplitMetricEnum>		SplitMetric(
				cl::desc("Metric for splitting:"),
				cl::values(clEnumValN(AverageRelVolume, "averV-metric", "Relative volume difference. Take average of MIV/MSV."),
						clEnumValN(AverageAbsVolume, "aveaV-metric", "Absolute volume difference. Take average of MIV/MSV."),
				clEnumValN(AverageHausdorff, "aveH-metric", "Hausdorff. Take average of MIV/MSV"),
				clEnumValN(SeparationVolume, "sepV-metric", "Volume in between resulting MIV/MSV"),
				clEnumValN(SharedOverlap, "shared-overlap-metric", "Vol difference between inbetween shape"),
				clEnumValEnd),cl::init(AverageRelVolume) );

enum ClusterMetricEnum
{
	SplitDist, CompleteLink, SingleLink
};

cl::opt<ClusterMetricEnum>		ClusterMetric(
				cl::desc("Metric for cluster packing distance:"),
				cl::values(clEnumValN(SplitDist, "split-cmetric", "Use splitting metric between MIV/MSV representations of clusters"),
						clEnumValN(CompleteLink, "complete-cmetric", "Use complete linkage value between cluster members"),
						clEnumValN(SingleLink, "single-cmetric", "Use single linkage value between cluster members"),
				clEnumValEnd),cl::init(SplitDist) );



extern cl::opt<bool> Verbose;

enum TopDownEnum
{
	OctantSplit, KSampleSplit
};

cl::opt<TopDownEnum>		TopDownSplit(
				cl::desc("Top down splitting method:"),
				cl::values(clEnumValN(OctantSplit, "octant-split", "Split based on informative octants"),
				clEnumValN(KSampleSplit, "ksample-split", "Spit using centers of random sample"),
				clEnumValEnd),cl::init(KSampleSplit) );

enum CenterFindEnum
{
	AveCenter,MinMaxCenter
};

cl::opt<CenterFindEnum>		CenterFind(
				cl::desc("Cluster center finding method:"),
				cl::values(clEnumValN(AveCenter, "ave-center", "Center defined by best average"),
				clEnumValN(MinMaxCenter, "minmax-center", "Center defined by minmax distance"),
				clEnumValEnd) ,cl::init(AveCenter));

enum PackingEnum
{
	IterativeMerge, FullMerge
};

cl::opt<PackingEnum>		PackingAlgorithm(
				cl::desc("Cluster packing method:"),
				cl::values(clEnumValN(IterativeMerge, "oitr-merge", "Greedy iterative merging"),
				clEnumValN(FullMerge, "ofull-merge", "Merge everything simultaneously"),
				clEnumValEnd) ,cl::init(FullMerge));


const unsigned GSSTree::MaxSplit = 8; //number of children in each node

extern cl::opt<unsigned> KCenters;
extern cl::opt<unsigned> KSampleMult;

GSSTree::GSSNode::~GSSNode()
{
	delete MSV;
	delete MIV;
}

void GSSTree::setBoundingBox(array<float, 6>& box)
{
	float dims[3] =
	{ 0, 0, 0 };
	for (unsigned i = 0; i < 3; i++)
	{
		min[i] = box[2 * i];
		dims[i] = box[2 * i + 1] - min[i];
	}

	//to make things easier, cubify using the longest dimension
	dim = max(dims[2], max(dims[0], dims[1]));
	//round all dimensions up to a power of 2; use fp
	dim = pow(2, ceil(log(dim) / log(2)));
	if (dim < 1)
		dim = 1;

	float div = dim;
	maxlevel = 0;
	while (div > maxres)
	{
		maxlevel++;
		div /= 2;
	}
}

// perform any transformations necessary to conver the input coordinate system
//into the gss tree coordinate system
void GSSTree::transformMol(const vector<MolSphere>& m, vector<MolSphere>& ret)
{
	ret.clear();
	ret.reserve(m.size());

	for (unsigned i = 0, n = m.size(); i < n; i++)
	{
		ret.push_back(
				MolSphere(m[i].x - min[0], m[i].y - min[1], m[i].z - min[2],
						m[i].r));
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
	if (KBestInsertion > 0)
	{
		root->findInsertionPoints(oct, kbest, KBestInsertion);

		assert(kbest.size() > 0);
		//find the one that creates the smallest change in the MIV/MSV of the leaf
		unsigned besti = 0;
		float smallestChange = HUGE_VAL;
		for (unsigned i = 0, n = kbest.size(); i < n; i++)
		{
			float change = kbest[i].leaf->combinedVolumeChange(oct, oct);
			if (change < smallestChange)
			{
				smallestChange = change;
				besti = i;
			}
		}
		assert(isfinite(smallestChange));
		leaf = kbest[besti].leaf;

		if (!KBestReshuffle)
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
double GSSTree::Partitioner::shannonEntropy(unsigned *patternCnts,
		unsigned total)
{
	double res = 0;
	double tot = total;
	for (unsigned i = 0; i < 256; i++)
	{
		double pi = patternCnts[i] / tot;
		if (pi > 0)
		{
			res += pi * log(pi);
		}
	}

	return -res;
}

//create additional octants at a higher level if need be
void GSSTree::Partitioner::updateOctantsForLevel(unsigned level,
		OctantCoords& octants, vector<bool>& done)
{
	//the number of octants grows exponentially with the level, for now
	//do the brute force approach, but eventually replace with a fraction ptr-based octtree
	//approach that traverse only the variable-presence octants
	if (octants.size() == 0)
	{
		octants.push_back(vector<unsigned> ());
		done.clear();
		done.resize(1, false);
	}

	if (level > octants[0].size())
	{
		vector<vector<unsigned> > newoctants;
		newoctants.reserve(octants.size() * 8);
		for (unsigned i = 0; i < 8; i++)
		{
			for (unsigned j = 0, n = octants.size(); j < n; j++)
			{
				newoctants.push_back(octants[j]);
				newoctants.back().push_back(i);
			}
		}
		octants.swap(newoctants);
		done.resize(octants.size(), false);
	}
}

// find an octant at the specified level that is "good"
// if splitMSV is true, look only at MSV patterns, otherwise MIV
// put result in octantcoord
bool GSSTree::Partitioner::findOctantAndSetUsed(unsigned level, bool splitMSV,
		vector<unsigned>& octantcoord)
{
	//some ideas for "best" octant
	//-largest variance in volume
	//-most populated bit patterns
	//-shannon index of bit patterns  <----
	//-other diversity indices of bit patterns

	OctantCoords& octants = MSVoctants;
	vector<bool>& done = MSVdoneOct;
	if (!splitMSV)
	{
		octants = MIVoctants;
		done = MIVdoneOct;
	}

	updateOctantsForLevel(level, octants, done);

	double bestval = 0;
	unsigned besto = 0;
	for (unsigned o = 0, no = octants.size(); o < no; o++)
	{
		if (done[o])
			continue;
		unsigned patternCnts[256] =
		{ 0, };
		unsigned total = 0;
		for (unsigned i = 0, n = tindex.size(); i < n; i++)
		{
			unsigned index = tindex[i];
			unsigned which = 0;
			if (splitMSV)
				which = partdata->getMSV(index)->getOctantPattern(octants[o],
						true);
			else
				which = partdata->getMIV(index)->getOctantPattern(octants[o],
						false);
			patternCnts[which]++;
			total++;
		}

		//various metrics are possible, probably just evaluate
		//shannon entropy and volume (faster?)
		double entropy = shannonEntropy(patternCnts, total);
		if (entropy > bestval)
		{
			bestval = entropy;
			besto = o;
		}
	}

	if (bestval == 0)
		return false; //all the same at this level

	octantcoord = octants[besto];
	done[besto] = true;
	return true;
}

//create sub partitions based on the octant bit pattern at the specified coord
void GSSTree::Partitioner::partitionOnOctant(
		const vector<unsigned>& octantcoord, bool splitMSV,
		vector<Partitioner>& parts)
{
	int pos[256];
	memset(pos, -1, sizeof(pos));

	parts.clear();
	parts.reserve(256);
	for (unsigned i = 0, n = tindex.size(); i < n; i++)
	{
		unsigned index = tindex[i];
		unsigned which = 0;
		if (splitMSV)
			which
					= partdata->getMSV(index)->getOctantPattern(octantcoord,
							true);
		else
			which = partdata->getMIV(index)->getOctantPattern(octantcoord,
					false);

		if (pos[which] < 0) //create new partition
		{
			pos[which] = parts.size();
			parts.push_back(Partitioner());
			parts.back().inheritFrom(*this);
		}

		parts[pos[which]].addSingle(*this, i);
	}
}

//merge from into this, must have the same data
//does no checking for duplicates
void GSSTree::Partitioner::add(const Partitioner& from)
{
	assert(partdata == from.partdata);
	for (unsigned i = 0, n = from.size(); i < n; i++)
		tindex.push_back(from.tindex[i]);
}

//add a single bit of partition data from from, updating structures as necessary
void GSSTree::Partitioner::addSingle(const Partitioner& from,
		unsigned fromindex)
{
	unsigned i = from.tindex[fromindex];
	tindex.push_back(i);

	if (MSV == NULL)
		MSV = partdata->getMSV(i)->clone();
	else
		MSV->unionWith(partdata->getMSV(i));

	if (MIV == NULL)
		MIV = partdata->getMIV(i)->clone();
	else
		MIV->intersect(partdata->getMIV(i));
}

void GSSTree::Partitioner::moveClusterInto(const Partitioner& from, GSSTree::Partitioner::Cluster& c)
{
	tindex.clear();
	for(unsigned i = 0, n = c.indices.size(); i < n; i++)
	{
		tindex.push_back(from.tindex[c.indices[i]]);
	}

	if(MSV) delete MSV;
	if(MIV) delete MIV;

	MIV = c.MIV;
	MSV = c.MSV;
	c.MIV = NULL;
	c.MSV = NULL;
	c.clear();
}

//put everything in rhs into this
//assumes octants are the same
void GSSTree::Partitioner::mergeWith(const Partitioner& rhs)
{
	assert(rhs.partdata == partdata);

	for(unsigned i = 0, n = rhs.tindex.size(); i < n; i++)
	{
		addSingle(rhs, i);
	}

	assert(MSVdoneOct.size() == rhs.MSVdoneOct.size());
	for(unsigned i = 0, n = MSVdoneOct.size(); i < n; i++)
	{
		MSVdoneOct[i] = MSVdoneOct[i] | rhs.MSVdoneOct[i];
	}
	assert(MIVdoneOct.size() == rhs.MIVdoneOct.size());
	for(unsigned i = 0, n = MIVdoneOct.size(); i < n; i++)
	{
		MIVdoneOct[i] = MIVdoneOct[i] | rhs.MIVdoneOct[i];
	}
}


float GSSTree::Partitioner::Cluster::distance(const GSSTree::Partitioner::Cluster& rhs, const GSSTree::Partitioner& part,const boost::multi_array<float, 2>& Dcache) const
{
	switch(ClusterMetric)
	{
	case SplitDist:
		return splitDist(MIV,MSV, rhs.MIV, rhs.MSV);
	case CompleteLink:
	{
		//this is the maximum of the minimum distances between cluster members
		float max = 0;
		for(unsigned i = 0, ni = indices.size(); i < ni; i++)
		{
			float min = HUGE_VAL;
			for(unsigned j = 0, nj = rhs.indices.size(); j < nj; j++)
			{
				float dist = 0;
				if(Dcache.size() > 0) //usecache
				{
					dist = Dcache[indices[i]][rhs.indices[j]];
				}
				else //have tor ecompute
				{
					unsigned l = part.tindex[indices[i]];
					unsigned r = part.tindex[rhs.indices[j]];

					dist = splitDist(part.partdata->getMIV(l), part.partdata->getMSV(l),
						part.partdata->getMIV(r), part.partdata->getMSV(r));
				}

				if(dist < min)
					min = dist;
			}
			if(min > max)
				max = min;
		}
		return max;
	}
	case SingleLink:
	{
		//the minimum distance overall
		float min = HUGE_VAL;
		for(unsigned i = 0, ni = indices.size(); i < ni; i++)
		{
			for(unsigned j = 0, nj = rhs.indices.size(); j < nj; j++)
			{
				float dist = 0;
				if(Dcache.size() > 0) //usecache
				{
					dist = Dcache[indices[i]][rhs.indices[j]];
				}
				else //have tor ecompute
				{
					unsigned l = part.tindex[indices[i]];
					unsigned r = part.tindex[rhs.indices[j]];
					dist = splitDist(part.partdata->getMIV(l), part.partdata->getMSV(l),
						part.partdata->getMIV(r), part.partdata->getMSV(r));
				}

				if(dist < min)
					min = dist;
			}
		}
		return min;
	}
	}
	return 0;
}

//distances between i and j
struct IntraClusterDist
{
	unsigned i;
	unsigned j;
	float dist;

	IntraClusterDist(): i(0), j(0), dist(HUGE_VAL) {}

	IntraClusterDist(unsigned I, unsigned J, float d): i(I), j(J), dist(d) {}

	bool operator<(const IntraClusterDist& rhs) const
	{
		return dist < rhs.dist;
	}
};

typedef bool (*ICDComp)(const IntraClusterDist& lhs, const IntraClusterDist& rhs);
static bool reverseDist(const IntraClusterDist& lhs, const IntraClusterDist& rhs)
{
	return lhs.dist > rhs.dist;
}

//split this partition into clusters that are within maxsplit size; turn these into
//smaller partitions
void GSSTree::Partitioner::packClusters(vector<Partitioner>& clusters)
{
	if(tindex.size() == 0)
		return;
	vector<Cluster> clusts; clusts.resize(tindex.size());
	for(unsigned i = 0, n = tindex.size(); i < n; i++)
	{
		unsigned index = tindex[i];
		clusts[i].set(i, partdata->getMIV(index), partdata->getMSV(index));
	}

	switch(PackingAlgorithm)
	{
	case FullMerge:
		{
			//combine everything as much as possible
			while(fullMergeClusters(clusts))
				;
		}
		break;
	case IterativeMerge:
		iterativeMergeClusters(clusts);
		break;
	}

	//now create partitions
	clusters.clear();
	clusters.reserve(clusts.size());
	for (unsigned i = 0, n = clusts.size(); i < n; i++)
	{
		clusters.push_back(Partitioner(partdata, maxlevel));
		clusters.back().moveClusterInto(*this, clusts[i]);
	}
}

//a priority q for distances, roll our own to support more efficient
//removal of clusters and re-use of distances array
class DistancePQ
{
	vector<IntraClusterDist>& Q;

public:
	DistancePQ(vector<IntraClusterDist>& q): Q(q)
	{
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

	bool empty() const { return Q.size() == 0; }

	const IntraClusterDist& top() const { return Q[0]; }

	void push(const IntraClusterDist& c)
	{
		Q.push_back(c);
		push_heap(Q.begin(), Q.end(), reverseDist);
	}

	void pop()
	{
		pop_heap(Q.begin(), Q.end(), reverseDist);
		Q.pop_back();
	}

	//remove all references to clusters i and j; this is linear
	void removeIJ(unsigned i, unsigned j)
	{
		unsigned c = 0;
		unsigned end = Q.size();
		while(c < end)
		{
			if(Q[c].i != i && Q[c].i != j && Q[c].j != i && Q[c].j != j)
			{
				c++;
			}
			else
			{
				end--;
				swap(Q[c],Q[end]);
			}
		}

		Q.erase(Q.begin()+end, Q.end());
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

};

//this will greedly merge clusters, unlike full merge, it does not try to merge
//everything at once, but greedily builds close clusters
//ensures that all clusters have at least two elements
void GSSTree::Partitioner::iterativeMergeClusters(vector<GSSTree::Partitioner::Cluster>& clusters)
{
	unsigned N = clusters.size();
	clusters.reserve(N*2);

	//compute all intra cluster distances
	vector<IntraClusterDist> distances; distances.reserve(N*N / 2);
	boost::multi_array<float, 2> darray(boost::extents[N][N]);

	for (unsigned i = 0; i < N; i++)
	{
		darray[i][i] = 0;
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusters[i].distance(clusters[j], *this);
			darray[i][j] = darray[j][i] = dist;
			distances.push_back(IntraClusterDist(i,j,dist));
		}
	}

	DistancePQ pQ(distances); //manipulates a reference of distances

	unsigned max = MaxSplit;
	while(!pQ.empty())
	{
		//get the smallest distance
		unsigned i = pQ.top().i;
		unsigned j = pQ.top().j;
		pQ.pop();

		if(clusters[i].isValid() && clusters[j].isValid())
		{
			if(clusters[i].size() + clusters[j].size() <= max)
			{
				unsigned c = clusters.size();
				//merge, create a new cluster
				clusters.push_back(Cluster());
				clusters.back().mergeInto(clusters[i], clusters[j]);

				pQ.removeIJ(i, j);
				//now compute the distance between this cluster and all remaining clusters
				//and add to priority Q
				for(unsigned d = 0; d < c; d++)
				{
					if(clusters[d].isValid())
					{
						float dist = clusters[d].distance(clusters.back(), *this, darray);
						pQ.push(IntraClusterDist(c,d,dist));
					}
				}
			}
			else if(AdaptivePacking && max > 2)
			{
				max--;
			}
		}
	}

	//remove empty clusters
	vector<Cluster> newclusters; newclusters.reserve(clusters.size());
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		if(clusters[i].isValid())
		{
			newclusters.push_back(Cluster());
			newclusters.back().moveInto(clusters[i]);
		}
	}

	swap(clusters,newclusters);
}

//each cluster has one or more indices into tindex as well as the MSV/MIV for the cluster
//merge everything once, return true if did any merges
bool GSSTree::Partitioner::fullMergeClusters(vector<GSSTree::Partitioner::Cluster >& clusters)
{
	unsigned N = clusters.size();
	if(N == 1)
		return false;

	vector<IntraClusterDist> distances; distances.reserve(N*N / 2);
	unsigned maxclust = 0;

	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusters[i].distance(clusters[j], *this);
			distances.push_back(IntraClusterDist(i,j,dist));
		}
	}
	sort(distances.begin(), distances.end());

	vector<Cluster> newclusters; newclusters.reserve(N/2);
	unsigned merged = 0;
	bool didmerge = false;
	for(unsigned d = 0, maxdists = distances.size(); d < maxdists; d++)
	{
		//merge if not already merged
		unsigned i = distances[d].i;
		unsigned j = distances[d].j;
		if (clusters[i].isValid() && clusters[j].isValid()
				&& clusters[i].size() + clusters[j].size() <= MaxSplit)
		{
			newclusters.push_back(Cluster());
			newclusters.back().mergeInto(clusters[i], clusters[j]);
			merged += 2;
			didmerge = true;
			if (newclusters.back().size() > maxclust)
				maxclust = newclusters.back().size();
		}

		if(merged >= N-1)
			break;
	}

	if(merged != N) //find the loner(s), keep it as a singleton (may also be too big)
	{
		for (unsigned i = 0; i < N; i++)
		{
			if(clusters[i].isValid())
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
		}
	}


	//this method may leave a singleton at then end, merge it
	if(!didmerge && newclusters.size() > 1 && newclusters.back().size() == 1)
	{
		//push the lone singleton into some other cluster
		double minval = HUGE_VAL;
		unsigned best = 0;
		for(unsigned i = 0, n = newclusters.size() - 1; i < n; i++)
		{
			float dist = newclusters.back().distance(newclusters[i], *this);
			if(dist < minval)
			{
				minval = dist;
				best = i;
			}
		}

		newclusters[best].addInto(newclusters.back());
		newclusters.pop_back();
	}

	swap(newclusters,clusters);

	return didmerge;
}

//use aglomerative clusters (complete linkage) to create k clusters of variable size
void GSSTree::Partitioner::kClusters(unsigned k,
		vector<Partitioner>& clusters)
{
	//compute pairwise distances between all members
	unsigned N = tindex.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);
	vector<vector<unsigned> > clusts;
	clusts.reserve(1 + size());
	for (unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const OctTree *imiv = partdata->getMIV(tindex[i]);
		const OctTree *imsv = partdata->getMSV(tindex[i]);
		for (unsigned j = 0; j < i; j++)
		{
			const OctTree *jmiv = partdata->getMIV(tindex[j]);
			const OctTree *jmsv = partdata->getMSV(tindex[j]);
			if (imiv == imsv && jmiv == jmsv) //leaf case
				distances[i][j] = distances[j][i] = leafDist(imiv, jmiv);
			else
				distances[i][j] = distances[j][i] = splitDist(imiv, imsv, jmiv,
						jmsv);
		}
		clusts.push_back(vector<unsigned> ());
		clusts.back().push_back(i);
	}

	//iteratively merge all clusters until we run into size limit
	while (clusts.size() > k)
	{
		//find two clusters with smallest complete linkage distance
		//(minimize the maximum distance)
		float mindist = HUGE_VAL;
		unsigned besti = 0;
		unsigned bestj = 0;
		for (unsigned i = 0, n = clusts.size(); i < n; i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				float maxdist = 0;
				//find the max distance between i and j
				for (unsigned I = 0, ni = clusts[i].size(); I < ni; I++)
				{
					for (unsigned J = 0, nj = clusts[j].size(); J < nj; J++)
					{
						float d = distances[clusts[i][I]][clusts[j][J]];
						if (d > maxdist)
						{
							maxdist = d;
						}
					}
				}
				if (maxdist < mindist)
				{
					mindist = maxdist;
					besti = i;
					bestj = j;
				}
			}
		}

		if (besti == bestj)
			break; //couldn't merge any

		//move besti and bestj to end of vector
		unsigned last = clusts.size() - 1;
		if (bestj == last)
		{
			swap(clusts[besti], clusts[last - 1]);
		}
		else
		{
			swap(clusts[besti], clusts[last]);
			swap(clusts[bestj], clusts[last - 1]);
		}

		//insert into second to last
		clusts[last - 1].insert(clusts[last - 1].end(), clusts[last].begin(),
				clusts[last].end());
		clusts.pop_back(); //remove
	}

	if (Verbose)
		cout << "  pack " << tindex.size() << " " << clusts.size() << "\n";

	//now create partitions
	clusters.clear();
	clusters.reserve(clusts.size());
	for (unsigned i = 0, n = clusts.size(); i < n; i++)
	{
		clusters.push_back(Partitioner(partdata, maxlevel));
		for (unsigned j = 0, m = clusts[i].size(); j < m; j++)
		{
			clusters.back().addSingle(*this, clusts[i][j]);
		}
	}
}

//cluster a set of partitions, merging each cluster into a single partition
//tries to produce a balanced set of clusters, where not cluster is
//(much?) larger than the largest incoming partition
void GSSTree::Partitioner::clusterPartitions(vector<Partitioner>& clusters)
{
	if (clusters.size() <= 1)
		return;
	unsigned max = 0;
	//find the largest cluster
	unsigned largestpos = 0;
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		if (clusters[i].size() > max)
		{
			max = clusters[i].size();
			largestpos = i;
		}
	}
	max = std::min(16U,max);

	//compute pairwise distances between mergeable clusters
	//(combined size <= max)
	boost::multi_array<float, 2> distances(
			boost::extents[clusters.size()][clusters.size()]);
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		distances[i][i] = HUGE_VAL;
		for (unsigned j = 0; j < i; j++)
		{
			if (clusters[i].size() + clusters[j].size() >= max)
				distances[i][j] = distances[j][i] = HUGE_VAL;
			else
			{
				distances[i][j] = distances[j][i] = splitDist(
						clusters[i].getMIV(), clusters[i].getMSV(),
						clusters[j].getMIV(), clusters[j].getMSV());
			}
		}
	}

	while(clusters.size() > 1)
	{
		//find closest partitions
		float min = HUGE_VAL;
		int posi = -1, posj = -1;
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			for(unsigned j = 0; j < i; j++)
			{
				if(distances[i][j] < min)
				{
					min = distances[i][j];
					posi = i;
					posj = j;
				}
			}
		}
		//if none, done
		if(posi < 0)
			break;

		//merge
		clusters[posi].mergeWith(clusters[posj]);
		clusters[posj].clear();
		//mark posj as unmergeable
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			distances[posj][i] = distances[i][posj] = HUGE_VAL;
		}
		//recompute distances from merged cluster, old j is now unmergeable
		for(unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			if(!isfinite(distances[posi][i]))
				; //stays unmergeable
			else if (clusters[posi].size() + clusters[i].size() >= max)
				distances[posi][i] = distances[i][posi] = HUGE_VAL;
			else
			{
				distances[posi][i] = distances[i][posi] = splitDist(
						clusters[i].getMIV(), clusters[i].getMSV(),
						clusters[posi].getMIV(), clusters[posi].getMSV());
			}
		}
	}
}

//create a leaf node from the contents of partitioner, does not check size
GSSTree::GSSLeafNode* GSSTree::leafFromPartition(Partitioner& partitioner,
		LeafPartitionData& leafdata)
{
	GSSLeafNode *newnode = new GSSLeafNode(octGen, dim, maxres);

	for (unsigned i = 0, n = partitioner.size(); i < n; i++)
	{
		newnode->trees.push_back(leafdata.getTree(partitioner.getDataIndex(i)));
		newnode->data.push_back(leafdata.getData(partitioner.getDataIndex(i)));
	}
	newnode->selfUpdate(); //update MIV/MSV

	return newnode;
}

//create an internal node from the contents of partitioner, does not check size
GSSTree::GSSInternalNode* GSSTree::nodeFromPartition(Partitioner& partitioner,
		NodePartitionData& nodedata)
{
	GSSInternalNode *newnode = new GSSInternalNode(octGen, dim, maxres);

	for (unsigned i = 0, n = partitioner.size(); i < n; i++)
	{
		newnode->children.push_back(
				nodedata.getNode(partitioner.getDataIndex(i)));
	}
	newnode->selfUpdate(); //update MIV/MSV

	return newnode;
}

//shouldnt top down partition, fall back to bottom up
bool GSSTree::Partitioner::unableToPartition() const
{
	if(level > maxlevel)
		return true;
	if(TopDownSplit == KSampleSplit)
	{
		if(size() < KCenters)
			return true;
	}

	return false;
}

//find the value that is most central in this partition
//this is assumed to be a small partition since we calculate all n^2 distances
void GSSTree::Partitioner::getCenter(const OctTree *& MIV, const OctTree *& MSV) const
{
	unsigned N = tindex.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);

	//compute distances
	for (unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const OctTree *imiv = partdata->getMIV(tindex[i]);
		const OctTree *imsv = partdata->getMSV(tindex[i]);
		for (unsigned j = 0; j < i; j++)
		{
			const OctTree *jmiv = partdata->getMIV(tindex[j]);
			const OctTree *jmsv = partdata->getMSV(tindex[j]);
			distances[i][j] = distances[j][i] = splitDist(imiv, imsv, jmiv,
						jmsv);
		}
	}

	//what is most central?  the row with the lowest average?
	//or the minimum maximum value?
	float bestave = HUGE_VAL;
	unsigned bestavei = 0;
	float minmaxval = HUGE_VAL;
	unsigned minmaxi = 0;
	for (unsigned i = 0; i < N; i++)
	{
		float ave = 0;
		float max = 0;
		for(unsigned j = 0; j < N; j++)
		{
			ave += distances[i][j];
			if(distances[i][j] > max)
				max = distances[i][j];
		}
		ave /= N;

		if(ave < bestave)
		{
			bestave = ave;
			bestavei = i;
		}
		if(max < minmaxval)
		{
			minmaxval = max;
			minmaxi = i;
		}
	}

	unsigned best;
	if(CenterFind == AveCenter)
		best = bestavei;
	else
		best = minmaxi;

	MIV = partdata->getMIV(tindex[best]);
	MSV = partdata->getMSV(tindex[best]);
}

//performs a top down partition, select a random sample, clusters it to
//get k clusters, and then uses the cluster centers to partition the whole data set
void GSSTree::Partitioner::topDownKSamplePartition(vector<Partitioner>& parts)
{
	//grab k*mult samples
	unsigned nsamples = KCenters*KSampleMult;
	assert(nsamples > 0);
	//assume reasonable input distribution, just slice out instead of random
	unsigned inc = size()/nsamples;
	if(inc == 0) inc = 1;

	Partitioner samples;
	samples.inheritFrom(*this);
	for(unsigned i = 0, n = size(); i < n; i += inc)
	{
		samples.addSingle(*this, i);
	}

	vector<Partitioner> clusters;
	samples.kClusters(KCenters, clusters);

	vector<const OctTree*> MSVcenters; MSVcenters.reserve(KCenters);
	vector<const OctTree*> MIVcenters; MIVcenters.reserve(KCenters);

	parts.resize(clusters.size());
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		MSVcenters.push_back(NULL);
		MIVcenters.push_back(NULL);

		clusters[i].getCenter(MIVcenters.back(), MSVcenters.back());
		parts[i].inheritFrom(*this);
	}

	unsigned numclusters = clusters.size();
	//now add each data iterm to the cluster whose center it is closest to
	for (unsigned i = 0, n = tindex.size(); i < n; i++)
	{
		unsigned index = tindex[i];

		//TODO: use triangle inequality to make this more efficient
		float mindist = HUGE_VAL;
		unsigned best = 0;
		for(unsigned j = 0; j < numclusters; j++)
		{
			float d = splitDist(MIVcenters[j],MSVcenters[j], partdata->getMIV(index), partdata->getMSV(index));
			if(d < mindist)
			{
				mindist = d;
				best = j;
			}
		}

		parts[best].addSingle(*this, i);
	}

	//identify any singletons and add them to an alternative cluster
	if (parts.size() > 1)
	{
		for (unsigned i = 0, n = parts.size(); i < n; i++)
		{
			if (parts[i].size() == 1)
			{
				unsigned index = parts[i].tindex.front();
				float mindist = HUGE_VAL;
				unsigned best = 0;
				for (unsigned j = 0; j < numclusters; j++)
				{
					if (j == i || parts[j].size() == 0)
						continue;
					float d = splitDist(MIVcenters[j],MSVcenters[j], partdata->getMIV(index), partdata->getMSV(index));
					if(d < mindist)
					{
						mindist = d;
						best = j;
					}
				}

				parts[best].addSingle(parts[i], 0);
				parts[i].clear();
			}
		}
	}

	if(Verbose)
	{
		cout << "KSample (" << size() << ") ";
		for(unsigned i = 0, n = parts.size(); i < n; i++)
		{
			cout << " " << parts[i].size();
		}
		cout << "\n";
	}
}

//performs a top down partition, this must be efficient O(n) and will
//create sub partitions for further partitioning while generating any
//small clusters that arise
void GSSTree::Partitioner::topDownOctantPartition(vector<Partitioner>& subparts)
{
	//identify octant to split on
	vector<unsigned> octantcoord;
	while (level <= maxlevel)
	{
		if (findOctantAndSetUsed(level, splitMSV, octantcoord))
			break;
		splitMSV = !splitMSV;
		if (findOctantAndSetUsed(level, splitMSV, octantcoord))
			break;
		level++;
	}

	//first partition based on the occupancy pattern of the specified octant
	vector<Partitioner> parts;
	partitionOnOctant(octantcoord, splitMSV, parts);

	if(PackPartitions)
		Partitioner::clusterPartitions(parts);

	//select which bins to recursively split and which to combine
	//into a single group to be clustered into leaves
	vector<unsigned> partitionMore;
	Partitioner mergedPartition(partdata, maxlevel);

	for (unsigned i = 0, n = parts.size(); i < n; i++)
	{
		if (parts[i].size() == 0)
			continue;
		if (parts[i].size() < LeafMerge)
		{
			mergedPartition.add(parts[i]);
			parts[i].clear();
		}
	}

	//create small partitions from small groups
	vector<Partitioner> clusters;
	mergedPartition.packClusters(clusters);

	//now fill out subparts with all the partitions
	subparts.clear(); subparts.reserve(parts.size() + clusters.size());
	for (unsigned i = 0, n = parts.size(); i < n; i++)
	{
		if(parts[i].size() > 0)
		{
			subparts.push_back(Partitioner());
			swap(subparts.back(),parts[i]);
		}
	}
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		subparts.push_back(Partitioner());
		swap(subparts.back(),clusters[i]);
	}
}

/*
 * Take all the trees indexed by tindex, segregate them
 */
void GSSTree::partitionLeaves(Partitioner& partitioner,
		LeafPartitionData& leafdata, vector<GSSNode*>& nodes)
{
	if (partitioner.size() == 0)
		return;

	if (partitioner.size() <= MaxSplit)
	{
		//create a node
		nodes.push_back(leafFromPartition(partitioner, leafdata));
		return;
	}


	//check for cases where we should just agglomerate
	if (partitioner.size() < LeafPack || partitioner.unableToPartition())
	{
		vector<Partitioner> clusters;
		partitioner.packClusters(clusters);
		for (unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(leafFromPartition(clusters[i], leafdata));
		}
	}
	else
	{
		//top down partition
		vector<Partitioner> parts;
		if(TopDownSplit == OctantSplit)
			partitioner.topDownOctantPartition(parts);
		else
			partitioner.topDownKSamplePartition(parts);
		//split on octants with large partitions
		for (unsigned i = 0, n = parts.size(); i < n; i++)
		{
			partitionLeaves(parts[i], leafdata, nodes);
		}
	}
}

//partition nodes to build another level
//very similar to partition leaves, but different data structures
//also, opportunaty for different algorithmic choices
void GSSTree::partitionNodes(Partitioner& partitioner,
		NodePartitionData& nodedata, vector<GSSNode*>& nodes)
{
	if (partitioner.size() == 0)
		return;
	if (partitioner.size() <= MaxSplit)
	{
		//create a node
		nodes.push_back(nodeFromPartition(partitioner, nodedata));
		return;
	}


	//check for cases where we should just agglomerate
	if (partitioner.size() < NodePack || partitioner.unableToPartition())
	{
		vector<Partitioner> clusters;
		partitioner.packClusters(clusters);
		for (unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			nodes.push_back(nodeFromPartition(clusters[i], nodedata));
		}
	}
	else
	{
		//top down partition
		vector<Partitioner> parts;
		if(TopDownSplit == OctantSplit)
			partitioner.topDownOctantPartition(parts);
		else
			partitioner.topDownKSamplePartition(parts);		//split on octants with large partitions
		for (unsigned i = 0, n = parts.size(); i < n; i++)
		{
			partitionNodes(parts[i], nodedata, nodes);
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
		for (unsigned i = 0, n = mols.size(); i < n; i++)
		{
			vector<MolSphere> mol;
			transformMol(mols[i], mol);

			trees.push_back(octGen.newOctTree(dim, maxres, mol));
			data.push_back(LeafData(mol));
		}

		LeafPartitionData leafdata(&trees, &data);
		Partitioner leafpartition(&leafdata, maxlevel);
		leafpartition.initFromData();
		vector<GSSNode*> nodes;
		vector<unsigned> octant;
		partitionLeaves(leafpartition, leafdata, nodes);
		trees.clear(); //now stored in nodes
		data.clear();

		while (nodes.size() > 1)
		{
			NodePartitionData nodedata(&nodes);
			Partitioner nodepartition(&nodedata, maxlevel);
			nodepartition.initFromData();
			vector<GSSNode*> nextlevel;
			partitionNodes(nodepartition, nodedata, nextlevel);
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
	if (Verbose)
		cout << "LeavesChecked " << leavesChecked << " / " << root->numLeaves()
				<< "\tFitChecks " << fitChecks << " " << elapsed << "\n";
	swap(res, data.spheres);

	if (ScanCheck)
	{
		//check
		t.restart();
		leavesChecked = 0;
		fitChecks = 0;
		dist = HUGE_VAL;
		root->scanNearest(tree, dist, data);
		elapsed = t.elapsed();
		if (Verbose)
			cout << "LeavesScanned " << leavesChecked << " / "
					<< root->numLeaves() << "\tFitChecks " << fitChecks << " "
					<< elapsed << "\n";

		if (data.spheres != res)
		{
			cout << "Find and Scan differ!\n";
		}
	}

	delete tree;
}

void GSSTree::tree_range_search(const OctTree* smallTree,
		const OctTree* bigTree, vector<vector<MolSphere> >& res)
{
	leavesChecked = 0;
	fitChecks = 0;
	Timer t;
	res.clear();
	vector<LeafData> data;
	root->findTweeners(smallTree, bigTree, data);
	double elapsed = t.elapsed();
	if (Verbose)
		cout << "LeavesChecked " << leavesChecked << " / " << root->numLeaves()
				<< "\tFitChecks " << fitChecks << " " << elapsed << "\t"
				<< data.size() << "\n";

	for (unsigned i = 0, n = data.size(); i < n; i++)
	{
		res.push_back(data[i].spheres);
	}

	if (ScanCheck)
	{
		//check
		t.restart();
		leavesChecked = 0;
		fitChecks = 0;
		data.clear();
		root->scanTweeners(smallTree, bigTree, data);
		elapsed = t.elapsed();
		if (Verbose)
			cout << "LeavesScanned " << leavesChecked << " / "
					<< root->numLeaves() << "\tFitChecks " << fitChecks << " "
					<< elapsed << "\t" << data.size() << "\n";

		if (data.size() != res.size()) //lame but easy check
		{
			cout << "Find and Scan differ!\n";
		}
	}
}

void GSSTree::dc_search(const vector<MolSphere>& little,
		const vector<MolSphere>& big, vector<vector<MolSphere> >& res)
{
	vector<MolSphere> littleMol, bigMol;
	transformMol(little, littleMol);
	transformMol(big, bigMol);

	OctTree *smallTree = octGen.newOctTree(dim, maxres, littleMol);
	OctTree *bigTree = octGen.newOctTree(dim, maxres, bigMol);

	tree_range_search(smallTree, bigTree, res);

	delete smallTree;
	delete bigTree;
}

void GSSTree::inex_search(const vector<MolSphere>& inc,
		const vector<MolSphere>& exc, vector<vector<MolSphere> >& res)
{
	vector<MolSphere> littleMol, bigMol;
	transformMol(inc, littleMol);
	transformMol(exc, bigMol);

	OctTree* smallTree = octGen.newOctTree(dim, maxres, littleMol);
	OctTree* bigTree = octGen.newOctTree(dim, maxres, bigMol);
	bigTree->invert();

	tree_range_search(smallTree, bigTree, res);

	delete smallTree;
	delete bigTree;
}

//output leaf data
void GSSTree::LeafData::write(ostream& out) const
{
	unsigned n = spheres.size();
	out.write((char*) &n, sizeof(n));
	for (unsigned i = 0; i < n; i++)
	{
		out.write((char*) &spheres[i], sizeof(MolSphere));
	}
}

//and input
void GSSTree::LeafData::read(istream& in)
{
	unsigned n = 0;
	in.read((char*) &n, sizeof(n));
	spheres.resize(n);
	for (unsigned i = 0; i < n; i++)
	{
		in.read((char*) &spheres[i], sizeof(MolSphere));
	}
}

//base class output
void GSSTree::GSSNode::write(ostream& out) const
{
	streampos ret = out.tellp();

	//write out kind of node
	const GSSInternalNode* intnode =
			dynamic_cast<const GSSInternalNode*> (this);
	bool isInternal = intnode != NULL;

	out.write((char*) &isInternal, sizeof(isInternal));
	out.write((char*) &res, sizeof(res));
	out.write((char*) &which, sizeof(which));
	//write out MIV/MSV
	MIV->write(out);
	MSV->write(out);
}

GSSTree::GSSNode* GSSTree::GSSNode::readCreate(const OctTreeFactory& octGen,
		istream& in, GSSInternalNode *parPtr)
{
	bool isInternal = false;
	in.read((char*) &isInternal, sizeof(isInternal));

	if (isInternal)
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
	in.read((char*) &res, sizeof(res));
	in.read((char*) &which, sizeof(which));
	MIV->read(in);
	MSV->read(in);
}

void GSSTree::GSSInternalNode::write(ostream& out) const
{
	GSSNode::write(out);

	unsigned n = children.size();
	out.write((char*) &n, sizeof(n));

	for (unsigned i = 0; i < n; i++)
	{
		children[i]->write(out);
	}
}

void GSSTree::GSSInternalNode::read(const OctTreeFactory& octGen, istream& in,
		GSSInternalNode *parPtr)
{
	GSSNode::read(in, parPtr);

	unsigned n = 0;
	in.read((char*) &n, sizeof(n));
	children.resize(n, NULL);
	for (unsigned i = 0; i < n; i++)
	{
		children[i] = GSSNode::readCreate(octGen, in, this);
	}
}

//dump leaf node
void GSSTree::GSSLeafNode::write(ostream& out) const
{
	GSSNode::write(out);

	unsigned n = trees.size();
	out.write((char*) &n, sizeof(n));
	for (unsigned i = 0; i < n; i++)
	{
		trees[i]->write(out);
	}

	for (unsigned i = 0; i < n; i++)
	{
		data[i].write(out);
	}
}

void GSSTree::GSSLeafNode::read(const OctTreeFactory& octGen, istream& in,
		GSSInternalNode *parPtr)
{
	GSSNode::read(in, parPtr);

	unsigned n = 0;
	in.read((char*) &n, sizeof(n));
	trees.resize(n, NULL);
	data.resize(n);

	for (unsigned i = 0; i < n; i++)
	{
		trees[i] = octGen.newOctTree();
		trees[i]->read(in);
	}

	for (unsigned i = 0; i < n; i++)
	{
		data[i].read(in);
	}
}

//dump tree in binary to be read into memory later
void GSSTree::write(ostream& out)
{
	out.write((char*) &maxres, sizeof(maxres));
	out.write((char*) min, sizeof(min));
	out.write((char*) &dim, sizeof(dim));
	octGen.write(out);
	root->write(out);
}

//read completely into memory
void GSSTree::read(istream& in)
{
	in.read((char*) &maxres, sizeof(maxres));
	in.read((char*) min, sizeof(min));
	in.read((char*) &dim, sizeof(dim));
	octGen.read(in);
	if (root)
		delete root;
	root = GSSNode::readCreate(octGen, in, NULL);
}

//evaluate the goodness of moving to from->to
//this is _expensive_ as we compute the change in volume
//for MIV/MSV of both leaves
float GSSTree::deltaFit(const GSSLeafNode *to, GSSLeafNode *from, unsigned t)
{
	//never remove the last tree from a node
	if (from->trees.size() <= 1)
		return 0;

	//compute MIV/MSV of from node without t
	OctTree* nMIV = octGen.newOctTree(dim, maxres);
	OctTree* nMSV = octGen.newOctTree(dim, maxres);
	nMIV->fill();

	for (unsigned i = 0, n = from->trees.size(); i < n; i++)
	{
		if (i != t)
		{
			nMIV->intersect(from->trees[i]);
			nMSV->unionWith(from->trees[i]);
		}
	}

	float fromMIV = nMIV->volume() - from->MIV->volume();
	float fromMSV = from->MSV->volume() - nMSV->volume();

	float fromCombined = fromMIV + fromMSV;
	float toCombined = to->combinedVolumeChange(from->trees[t], from->trees[t]);

	delete nMIV;
	delete nMSV;

	//only favorable if we tighten from more than we expand to
	return fromCombined - toCombined;
}

//insert data into leaf node, if there isn't enough room, split
//call update on parents if need-be
//takes ownership of the tree memory
void GSSTree::GSSLeafNode::insert(GSSTree& gTree, OctTree* tree,
		const LeafData& m, vector<LeafDistPair>& kbest)
{
	if (data.size() < MaxSplit)
	{
		bool needUpdate = false;

		assert(data.size() == trees.size());
		//easy case, just add
		trees.push_back(tree);
		data.push_back(LeafData(m));

		//update MIV/MSV
		needUpdate |= MIV->intersect(tree);
		needUpdate |= MSV->unionWith(tree);

		if (needUpdate && parent != NULL)
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
		GSSLeafNode *newnode = new GSSLeafNode(gTree.octGen, gTree.dim,
				gTree.maxres);

		//put split2 in new node
		newnode->MIV->fill();
		newnode->data.resize(split2.size());
		for (unsigned i = 0, n = split2.size(); i < n; i++)
		{
			unsigned index = split2[i];
			newnode->trees.push_back(trees[index]);
			swap(newnode->data[i], data[index]);

			//update MIV/MSV
			newnode->MSV->unionWith(trees[index]);
			newnode->MIV->intersect(trees[index]);
		}
		//now mogirfy to split1
		vector<OctTree*> tmptrees;
		tmptrees.reserve(MaxSplit);
		vector<LeafData> tmpdata(split1.size());

		MIV->fill();
		MSV->clear();
		for (unsigned i = 0, n = split1.size(); i < n; i++)
		{
			unsigned index = split1[i];
			tmptrees.push_back(trees[index]);
			swap(tmpdata[i], data[index]);

			MIV->intersect(trees[index]);
			MSV->unionWith(trees[index]);
		}

		swap(tmpdata, data);
		swap(tmptrees, trees);

		if (parent != NULL)
			parent->update(gTree, which, newnode);
		else if (newnode != NULL) //need new root
		{
			gTree.createRoot(this, newnode);
		}

		//if the kbest leaves were passed in, refine all the contents
		//of these leaves be evaluating whether or not each individual
		//tree would be more specific to one of the split nodes and moving it
		//if so - this is not cheap
		int cnt = ReshuffleLimit;
		if (kbest.size() > 0 && newnode != NULL)
		{
			for (unsigned i = 0, n = kbest.size(); i < n; i++)
			{
				if (cnt <= 0)
					break;
				GSSLeafNode *altleaf = kbest[i].leaf;
				if (altleaf == this)
					continue;

				//don't resplit, just stop if full
				if (trees.size() >= MaxSplit && newnode->trees.size()
						>= MaxSplit)
					break;

				bool changed1 = false, changed2 = false;
				for (int i = 0; i < (int) altleaf->trees.size(); i++)
				{
					float deltaFit1 = gTree.deltaFit(this, altleaf, i);
					float deltaFit2 = gTree.deltaFit(newnode, altleaf, i);

					//must be positive to be better fit
					if (deltaFit1 > 0 || deltaFit2 > 0)
					{
						if (deltaFit1 > deltaFit2 && trees.size() < MaxSplit)
						{
							//move to this node
							moveTreeFrom(altleaf, i);
							changed1 = true;
							i--; //has been deleted
							cnt--;
						}
						else if (newnode->trees.size() < MaxSplit)
						{
							//move to new node
							newnode->moveTreeFrom(altleaf, i);
							changed2 = true;
							i--;
							cnt--;
						}
					}
				}

				if (changed1 && parent != NULL) //update up the tree
					parent->fullUpdate();
				if (changed2 && newnode->parent != NULL)
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
	if (parent != NULL)
		parent->fullUpdate();
}

//update our MIV/MSV from scratch (necessary when removing)
//do not recurse upwards
void GSSTree::GSSLeafNode::selfUpdate()
{
	MIV->fill();
	MSV->clear();

	for (unsigned i = 0, n = trees.size(); i < n; i++)
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

	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		MIV->intersect(children[i]->MIV);
		MSV->unionWith(children[i]->MSV);
	}
}

//recurse UP the tree updating MIV/MSV and splitting as necessary
void GSSTree::GSSInternalNode::update(GSSTree& tree, unsigned whichChild,
		GSSNode *newnode)
{
	//first update MIV/MSV
	assert(whichChild < children.size());

	bool needUpdate = false;
	needUpdate |= MIV->intersect(children[whichChild]->MIV);
	needUpdate |= MSV->unionWith(children[whichChild]->MSV);

	if (newnode != NULL)
	{
		//add node to this level
		if (children.size() < MaxSplit)
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

			for (unsigned i = 0, n = children.size(); i < n; i++)
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
			for (unsigned i = 0, n = split1.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split1[i]];
				addChild(child);
				MIV->intersect(child->MIV);
				MSV->unionWith(child->MSV);
			}

			//same thing for the new node and split2
			GSSInternalNode *newinode = new GSSInternalNode(tree.octGen,
					tree.dim, tree.maxres);
			newinode->MIV->fill();
			newinode->MSV->clear();
			for (unsigned i = 0, n = split2.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split2[i]];
				newinode->addChild(child);
				newinode->MIV->intersect(child->MIV);
				newinode->MSV->unionWith(child->MSV);
			}
			newnode = newinode;
		}
	}

	if (parent != NULL)
	{
		parent->update(tree, which, newnode);
	}
	else if (newnode != NULL && parent == NULL) //need new root
	{
		tree.createRoot(this, newnode);
	}

}

//create a new root using left and right, which are assumed to
//be splits of the current root
void GSSTree::createRoot(GSSNode *left, GSSNode *right)
{
	assert(left == root);
	GSSInternalNode* newroot = new GSSInternalNode(octGen, dim, maxres);
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
bool GSSTree::fitsInbetween(const OctTree *MIV, const OctTree *MSV,
		const OctTree *min, const OctTree *max)
{
	fitChecks++;
	//TODO: more efficient
	//the MSV must completely enclose min
	if (!min->containedIn(MSV))
		return false;
	//MIV must be completely enclosed by max
	if (!MIV->containedIn(max))
		return false;

	return true;
}

//return a distance to a single leaf, should be compatible with search Dist
float GSSTree::leafDist(const OctTree* obj, const OctTree *leaf)
{
	switch(SplitMetric)
	{
	case AverageRelVolume:
	case SeparationVolume:
	case SharedOverlap:
		return obj->relativeVolumeDistance(leaf);
	case AverageAbsVolume:
		return obj->absoluteVolumeDistance(leaf);
		break;
	case AverageHausdorff:
		return max(obj->hausdorffDistance(leaf), leaf->hausdorffDistance(obj));

		break;
	}
	abort();
}

//return a "distance" between obj and MIV/MSV; the lower the distance
//the higher the more likely a node should be searched
//min and max should bookend the ultimate leaf distances
float GSSTree::searchDist(const OctTree* obj, const OctTree *MIV,
		const OctTree *MSV, float& min, float& max)
{
	min = 1 - obj->intersectVolume(MSV)/obj->unionVolume(MIV); //percent of shape already covered
	max = 1 - obj->intersectVolume(MIV)/obj->unionVolume(MSV); //difference if MIV/MSV after merge

	return min + max;
}

/*
 * return a "distance" between approximations
 * for now, just the sum of the volume different of the MIVs and of the MSVs
 */
float GSSTree::splitDist(const OctTree* leftMIV, const OctTree* leftMSV,
		const OctTree* rightMIV, const OctTree* rightMSV)
{
	//handle leaves differently
	if(leftMIV == leftMSV && rightMIV == rightMSV)
	{
		return leafDist(leftMIV, rightMIV);
	}

	switch(SplitMetric)
	{
	case AverageRelVolume:
		return leftMIV->relativeVolumeDistance(rightMIV) + leftMSV->relativeVolumeDistance(rightMSV);
	case AverageAbsVolume:
		return leftMIV->absoluteVolumeDistance(rightMIV) + leftMSV->absoluteVolumeDistance(rightMSV);
	case SeparationVolume:
		return leftMIV->inbetweenVolume(leftMSV, rightMIV, rightMSV);
		break;
	case SharedOverlap:
		return leftMIV->percentOverlapVolume(leftMSV, rightMIV, rightMSV);
	case AverageHausdorff:
		float d1 = max(leftMIV->hausdorffDistance(rightMIV),
						rightMIV->hausdorffDistance(leftMIV));
				float d2 = max(leftMSV->hausdorffDistance(rightMSV),
						rightMSV->hausdorffDistance(leftMSV));

				return d1 + d2;
		break;
	}
	abort();
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

	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = i + 1; j < N; j++)
		{
			float d = splitDist(MIV[i], MSV[i], MIV[j], MSV[j]);
			distances[i][j] = d;
			distances[j][i] = d;
			if (d > dist)
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

	for (unsigned i = 0; i < N; i++)
	{
		if (i != besti && i != bestj)
		{
			//choose based on distance to cumulative MIV/MSV
			//unless we've already filled one split
			if (s1.size() > N / 2 + 1)
			{
				s2.push_back(i);
			}
			else if (s2.size() > N / 2 + 1)
			{
				s1.push_back(i);
			}
			else
			{
				float di = splitDist(iMIV, iMSV, MIV[i], MSV[i]);
				float dj = splitDist(jMIV, jMSV, MIV[i], MSV[i]);

				if (di < dj)
				{
					s1.push_back(i);
					iMIV->intersect(MIV[i]);
					iMSV->unionWith(MSV[i]);
				}
				else if (dj < di)
				{
					s2.push_back(i);
					jMIV->intersect(MIV[i]);
					jMSV->unionWith(MSV[i]);
				}
				else //equal, seems unlikely
				{
					if (s1.size() < s2.size())
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
float GSSTree::GSSNode::combinedVolumeChange(const OctTree *miv,
		const OctTree *msv) const
{
	float deltamiv = MIV->volume() - MIV->intersectVolume(miv);
	float deltamsv = MSV->unionVolume(msv) - MSV->volume();
	return deltamiv + deltamsv;
}

//examine every object in the leaf and see if it has a better distance than distance,
//if so, store data
void GSSTree::GSSLeafNode::scanNearest(const OctTree* tree, float& distance,
		LeafData& d)
{
	findNearest(tree, distance, d); //this is just the same
}

//find the object with the best distance in this node, if it's better than
//the passed distance, update
void GSSTree::GSSLeafNode::findNearest(const OctTree* tree, float& distance,
		LeafData& d)
{
	leavesChecked++;
	for (unsigned i = 0, n = trees.size(); i < n; i++)
	{
		float dist = leafDist(tree, trees[i]);
		if (dist < distance)
		{
			distance = dist;
			d = data[i];
		}
	}
}

//identify all objects that are exactly in between min and max
void GSSTree::GSSLeafNode::scanTweeners(const OctTree* min, const OctTree* max,
		vector<LeafData>& res)
{
	findTweeners(min, max, res);
}

static unsigned foundLeafCnt = 0;
//identify all objects that are exactly in between min and max
void GSSTree::GSSLeafNode::findTweeners(const OctTree* min, const OctTree* max,
		vector<LeafData>& res)
{
	leavesChecked++;
	bool found = false;
	for (unsigned i = 0, n = trees.size(); i < n; i++)
	{
		if (fitsInbetween(trees[i], trees[i], min, max))
		{
			found = true;
			res.push_back(data[i]);
		}
	}

	if (0 && found)
	{
		//output stuff
		cout << "FOUNDLEAF " << foundLeafCnt << "\n";
		stringstream name;
		name << "LEAF" << foundLeafCnt << ".xyz";
		ofstream out(name.str().c_str());

		for (unsigned i = 0, n = trees.size(); i < n; i++)
		{
			bool good = fitsInbetween(trees[i], trees[i], min, max);
			cout << i << " " << good << " " << splitDist(trees[i], trees[i],
					min, max) << " |";
			for (unsigned j = 0; j < n; j++)
			{
				cout << " "
						<< splitDist(trees[j], trees[j], trees[i], trees[i]);
			}
			cout << "\n";

			//dump mol
			out << data[i].spheres.size();
			out << "\nLEAF" << foundLeafCnt << " " << i << "\n";
			for (unsigned j = 0, m = data[i].spheres.size(); j < m; j++)
			{
				out << "C " << data[i].spheres[j].x << " "
						<< data[i].spheres[j].y << " " << data[i].spheres[j].z
						<< "\n";
			}
		}
		cout << "----\n";
		foundLeafCnt++;

	}
}

//if this leaf is a more appropriate place for tree, return self
void GSSTree::GSSLeafNode::findInsertionPoint(const OctTree* tree,
		float& distance, GSSLeafNode*& leaf)
{
	if (trees.size() == 0)
	{
		distance = 0;
		leaf = this;
	}
	else
	{
		float min, max;
		float dist = searchDist(tree, MIV, MSV, min, max);
		if (dist < distance)
		{
			distance = dist;
			leaf = this;
		}
	}
}

//add this leaf if it's distance is better and truncate kbest if necessary
void GSSTree::GSSLeafNode::findInsertionPoints(const OctTree* tree,
		vector<LeafDistPair>& kbest, unsigned k)
{
	float dist = 0;
	if (trees.size() > 0)
	{
		float min, max;
		dist = searchDist(tree, MIV, MSV, min, max);
	}

	if (kbest.size() == 0 || kbest.back().distance > dist)
	{
		LeafDistPair pair(this, dist);
		kbest.insert(lower_bound(kbest.begin(), kbest.end(), pair), pair);
		if (kbest.size() > k) //truncate
			kbest.resize(k);
	}

}

//move the t'th tree of from into this, removing from from
//updates the local MIV/MSV, but does not recurse up the tree
void GSSTree::GSSLeafNode::moveTreeFrom(GSSLeafNode* from, unsigned t)
{
	assert(t < from->trees.size());

	swap(from->trees[t], from->trees.back());
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

	ScoreIndex()
	{
	}
	ScoreIndex(unsigned i, float s, float m) :
		index(i), score(s), min(m)
	{
	}

	bool operator<(const ScoreIndex& si) const
	{
		return score < si.score;
	}
};

//identify all objects that are exactly in between min and max
//brute force scan for debugging
void GSSTree::GSSInternalNode::scanTweeners(const OctTree* min,
		const OctTree* max, vector<LeafData>& res)
{
	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		children[i]->scanTweeners(min, max, res);
	}
}

//identify all objects that are exactly in between min and max
//filter out children that can't possibly match
void GSSTree::GSSInternalNode::findTweeners(const OctTree* min,
		const OctTree* max, vector<LeafData>& res)
{
	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		if (fitsInbetween(children[i]->MIV, children[i]->MSV, min, max))
		{
			children[i]->findTweeners(min, max, res);
		}
	}
}

//scan - ignore any bounding on the search
void GSSTree::GSSInternalNode::scanNearest(const OctTree* tree,
		float& distance, LeafData& data)
{
	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		children[i]->scanNearest(tree, distance, data);
	}
}

//explore children to find closest value
void GSSTree::GSSInternalNode::findNearest(const OctTree* tree,
		float& distance, LeafData& data)
{
	vector<ScoreIndex> SIs;
	SIs.reserve(children.size());

	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		float min = 0, max = 0;
		float score = searchDist(tree, children[i]->MIV, children[i]->MSV, min,
				max);
		if (min < distance) //there's hope
		{
			SIs.push_back(ScoreIndex(i, score, min));
		}
	}

	//explore in priority order
	sort(SIs.begin(), SIs.end());

	for (unsigned i = 0, n = SIs.size(); i < n; i++)
	{
		if (SIs[i].min < distance)
		{
			children[SIs[i].index]->findNearest(tree, distance, data);
		}
	}
}

//look for a good place for the data indexed by tree
void GSSTree::GSSInternalNode::findInsertionPoint(const OctTree* tree,
		float& distance, GSSLeafNode*& leaf)
{
	float min, max;
	float dist = searchDist(tree, MIV, MSV, min, max);
	if (dist < distance)
	{
		//if this tree might contain something better, look at all children -> really should sort this
		for (unsigned i = 0, n = children.size(); i < n; i++)
		{
			children[i]->findInsertionPoint(tree, distance, leaf);
		}
	}
}

void GSSTree::GSSInternalNode::findInsertionPoints(const OctTree* tree,
		vector<LeafDistPair>& kbest, unsigned k)
{
	float min, max;
	float dist = searchDist(tree, MIV, MSV, min, max);
	float bestdist = HUGE_VAL;
	if (kbest.size() == k && k > 0)
	{
		bestdist = kbest.back().distance;
	}
	if (dist < bestdist)
	{
		//if this tree might contain something better, look at all children -> really should sort this
		for (unsigned i = 0, n = children.size(); i < n; i++)
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
	if (sz > nodes.max)
		nodes.max = sz;
	if (sz < nodes.min)
		nodes.min = sz;
	nodes.total += sz;
	nodes.cnt++;
	if (sz == 1)
		nodes.singletonCnt++;
	unsigned ret = 0;
	for (unsigned i = 0, n = children.size(); i < n; i++)
	{
		ret = max(ret, children[i]->getStats(leaves, nodes));
	}
	return ret + 1;
}

unsigned GSSTree::GSSLeafNode::getStats(Stats& leaves, Stats& nodes) const
{
	unsigned sz = trees.size();
	if (sz > leaves.max)
		leaves.max = sz;
	if (sz < leaves.min)
		leaves.min = sz;
	leaves.total += sz;
	leaves.cnt++;
	if (sz == 1)
		leaves.singletonCnt++;
	return 0;
}

//print out aggregrate statistics for density of nodes/leaves
void GSSTree::printStats() const
{
	Stats leaves, nodes;
	unsigned depth = root->getStats(leaves, nodes);

	printf(
			"%dDepth %fAve | Leaves (%d): %fAve %dMin %dMax %dSingle| Nodes (%d): %fAve %dMin %dMax %dSingle\n",
			depth,
			(leaves.total + nodes.total) / (double) (leaves.cnt + nodes.cnt),
			leaves.cnt, leaves.total / (double) leaves.cnt, leaves.min,
			leaves.max, leaves.singletonCnt, nodes.cnt,
			nodes.total / (double) nodes.cnt, nodes.min, nodes.max,
			nodes.singletonCnt);
}

