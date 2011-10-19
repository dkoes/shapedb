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
//initialize to fully represent the data
TopDownPartitioner* KSamplePartitioner::create(const DataViewer* dv) const
{
	TopDownPartitioner* ret = new KSamplePartitioner(dv, kcenters, ksamples);
	ret->initFromData();
	return ret;
}

//initialize using a slice of the data from ind
TopDownPartitioner* KSamplePartitioner::create(const DataViewer* dv, vector<unsigned>& ind) const
{
	KSamplePartitioner* ret = new KSamplePartitioner(dv, kcenters, ksamples);
	swap(ret->indices, ind);
	return ret;
}



//find the SINGLE value that is most central in cluster
//this is assumed to be a small partition since we calculate all n^2 distances
void KSamplePartitioner::getCenter(const vector<unsigned>& cluster, const MappableOctTree *& MIV, const MappableOctTree *& MSV) const
{
	unsigned N = cluster.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);

	//compute distances
	for (unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const MappableOctTree *imiv = data->getMIV(cluster[i]);
		const MappableOctTree *imsv = data->getMSV(cluster[i]);
		for (unsigned j = 0; j < i; j++)
		{
			const MappableOctTree *jmiv = data->getMIV(cluster[j]);
			const MappableOctTree *jmsv = data->getMSV(cluster[j]);
			distances[i][j] = distances[j][i] = shapeDistance(imiv, imsv, jmiv,
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
	if(centerFind == AveCenter)
		best = bestavei;
	else
		best = minmaxi;

	MIV = data->getMIV(cluster[best]);
	MSV = data->getMSV(cluster[best]);
}


void KSamplePartitioner::partition(vector<TopDownPartitioner*>& parts)
{
	//grab k*mult samples
	unsigned nsamples = kcenters*ksamples;
	assert(nsamples > 0);
	//assume reasonable input distribution, just slice out instead of random
	unsigned inc = indices.size()/nsamples;
	if(inc == 0) inc = 1;

	vector<unsigned> sampleIndices; sampleIndices.reserve(nsamples+1);
	for(unsigned i = 0, n = indices.size(); i < n; i += inc)
	{
		sampleIndices.push_back(indices[i]);
	}

	//cluster samples
	vector< vector<unsigned> > clusters;
	kCluster(sampleIndices, clusters);

	//compute cluster "centers": MIV/MSV
	unsigned numclusters = clusters.size();

	vector<const MappableOctTree*> MSVcenters(numclusters, NULL);
	vector<const MappableOctTree*> MIVcenters(numclusters, NULL);
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		getCenter(clusters[i], MIVcenters[i], MSVcenters[i]);
	}

	vector< vector<unsigned> > partitions;
	partitions.resize(numclusters);

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

}


//use aglomerative clusters (complete linkage) to create k clusters of variable size
void KSamplePartitioner::kCluster(const vector<unsigned>& indices, vector< vector<unsigned> >& clusters)
{
	//compute pairwise distances between all members
	unsigned N = indices.size();
	boost::multi_array<float, 2> distances(boost::extents[N][N]);
	//compute distances indicesed by position in indices
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
		clusters.push_back(vector<unsigned>());
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

	//clusters is now populated with indicies into indices
	//replace with actual indices
	for(unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		for(unsigned j = 0, m = clusters[i].size(); j < m; j++)
		{
			clusters[i][j] = indices[clusters[i][j]];
		}
	}

}
