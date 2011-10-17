/*
 * KSamplePartitioner.cpp
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 */

#include "KSamplePartitioner.h"
#include "ShapeDistance.h"

#include <boost/multi_array.hpp>
using namespace boost;

//create an instances of ksamplepartitioner with all of this's settings and the specified data
TopDownPartitioner* KSamplePartitioner::create(const DataViewer* dv, const vector<file_index>& ind) const
{
	KSamplePartitioner* ret = new KSamplePartitioner(kcenters, ksamples);
	ret->data = dv;
	ret->indices = ind;
	return ret;
}

//take the intersection of everything represented by ind
MappableOctTree* KSamplePartitioner::computeMIV(const vector<file_index>& ind) const
{
	const MappableOctTree* trees[ind.size()];
	for(unsigned i = 0, n = ind.size(); i < n; i++)
	{
		trees[i] = data->getMIV(ind[i]);
	}
	return MappableOctTree::createFromIntersection(ind.size(), trees);
}

//take the union of everything represented by ind
MappableOctTree* KSamplePartitioner::computeMSV(const vector<file_index>& ind) const
{
	const MappableOctTree* trees[ind.size()];
	for(unsigned i = 0, n = ind.size(); i < n; i++)
	{
		trees[i] = data->getMSV(ind[i]);
	}
	return MappableOctTree::createFromUnion(ind.size(), trees);
}

void KSamplePartitioner::partition(vector<TopDownPartitioner*>& parts)
{
	//grab k*mult samples
	unsigned nsamples = kcenters*ksamples;
	assert(nsamples > 0);
	//assume reasonable input distribution, just slice out instead of random
	unsigned inc = indices.size()/nsamples;
	if(inc == 0) inc = 1;

	vector<file_index> sampleIndices; sampleIndices.reserve(nsamples+1);
	for(unsigned i = 0, n = indices.size(); i < n; i += inc)
	{
		sampleIndices.push_back(indices[i]);
	}

	//cluster samples
	vector< vector<file_index> > clusters;
	kCluster(sampleIndices, clusters);

	//compute cluster "centers": MIV/MSV
	vector<const MappableOctTree*> MSVcenters; MSVcenters.reserve(kcenters);
	vector<const MappableOctTree*> MIVcenters; MIVcenters.reserve(kcenters);

	vector< vector<file_index> > partitions;
	partitions.resize(clusters.size());
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		MSVcenters.push_back(computeMIV(clusters[i]));
		MIVcenters.push_back(computeMSV(clusters[i]));
	}

	unsigned numclusters = clusters.size();
	//now add each data iterm to the cluster whose center it is closest to
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];

		//TODO: use triangle inequality to make this more efficient
		float mindist = HUGE_VAL;
		unsigned best = 0;
		for(unsigned j = 0; j < numclusters; j++)
		{
			float d = shapeDistance(MIVcenters[j],MSVcenters[j], data->getMIV(index), data->getMSV(index));
			if(d < mindist)
			{
				mindist = d;
				best = j;
			}
		}

		partitions[best].push_back(index);
	}

	//identify any singletons and add them to an alternative cluster
	if (partitions.size() > 1)
	{
		for (unsigned i = 0, n = partitions.size(); i < n; i++)
		{
			if (partitions[i].size() == 1)
			{
				unsigned index = partitions[i].front();
				float mindist = HUGE_VAL;
				unsigned best = 0;
				for (unsigned j = 0; j < numclusters; j++)
				{
					if (j == i || partitions[j].size() == 0)
						continue;
					float d = shapeDistance(MIVcenters[j],MSVcenters[j], data->getMIV(index), data->getMSV(index));
					if(d < mindist)
					{
						mindist = d;
						best = j;
					}
				}

				partitions[best].push_back(index);
				partitions[i].clear();
			}
		}
	}

	//create new partitions from indices
	for (unsigned i = 0, n = partitions.size(); i < n; i++)
	{
		if(partitions[i].size() > 0)
		{
			parts.push_back(create(data, partitions[i]));
		}
	}

	//clear memory
	for(unsigned i = 0; i < numclusters; i++)
	{
		delete MIVcenters[i];
		delete MSVcenters[i];
	}
}


//use aglomerative clusters (complete linkage) to create k clusters of variable size
void KSamplePartitioner::kCluster(const vector<file_index>& indices, vector< vector<file_index> >& clusters)
{
	//compute pairwise distances between all members
	unsigned N = indices.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);
	clusters.clear();
	clusters.reserve(1 + N);
	for (unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const MappableOctTree *imiv = data->getMIV(indices[i]);
		const MappableOctTree *imsv = data->getMSV(indices[i]);
		for (unsigned j = 0; j < i; j++)
		{
			const MappableOctTree *jmiv = data->getMIV(indices[j]);
			const MappableOctTree *jmsv = data->getMSV(indices[j]);
			distances[i][j] = distances[j][i] = shapeDistance(imiv, imsv, jmiv, jmsv);
		}
		clusters.push_back(vector<file_index> ());
		clusters.back().push_back(i);
	}

	//iteratively merge all clusters until we run into size limit
	while (clusters.size() > kcenters)
	{
		//find two clusters with smallest complete linkage distance
		//(minimize the maximum distance)
		float mindist = HUGE_VAL;
		unsigned besti = 0;
		unsigned bestj = 0;
		for (unsigned i = 0, n = clusters.size(); i < n; i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				float maxdist = 0;
				//find the max distance between i and j
				for (unsigned I = 0, ni = clusters[i].size(); I < ni; I++)
				{
					for (unsigned J = 0, nj = clusters[j].size(); J < nj; J++)
					{
						float d = distances[clusters[i][I]][clusters[j][J]];
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
		unsigned last = clusters.size() - 1;
		swap(clusters[besti], clusters[last]);
		swap(clusters[bestj], clusters[last - 1]);

		//insert into second to last
		clusters[last - 1].insert(clusters[last - 1].end(), clusters[last].begin(),
				clusters[last].end());
		clusters.pop_back(); //remove
	}

}
