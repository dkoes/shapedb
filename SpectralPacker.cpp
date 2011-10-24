/*
 * SpectralPacker.cpp
 *
 *  Created on: Oct 24, 2011
 *      Author: dkoes
 */

#include "SpectralPacker.h"
#include "ShapeDistance.h"
#include <boost/shared_ptr.hpp>

//create a similarity matrix from graph
//compute log(n) nearest neighbors to determine a good sigma
//or alternatively compute an actual knn graph (todo - a little tricky
//since need to keep the graph connected)
void SpectralPacker::transformDistancesToSimilarity(const multi_array<double, 2>& distances, MatrixXd& graph) const
{
	unsigned N = distances.size(); //must be square
	graph.resize(N,N);
	unsigned k = max(ceil(log2(N)),1.0);

	vector<double> kthdist(N, 0);
	double sum = 0;
	//compute the k nearest
	vector<double> row(N);
	for(unsigned i = 0; i < N; i++)
	{
		//get distances
		for(unsigned j = 0; j < N; j++)
		{
			if(j == i)
				row[j] = HUGE_VAL;
			else
				row[j] = distances[i][j];
		}
		sort(row.begin(),row.end());
		kthdist[i] = row[k-1];
		sum += row[k-1];
	}

	double sigma = sum/N;
	double denom = 2*sigma*sigma;
	//now compute Gaussian similarity
	for(unsigned i = 0; i < N; i++)
	{
		for(unsigned j = 0; j < N; j++)
		{
			if(i == j)
				graph(i,j) = 0;
			else
			{
				double d = distances[i][j];
				d *= d; //sqr
				graph(i,j) = exp(-d/denom);
			}
		}
	}
}

//The unnormalized Laplacian is D - W
void SpectralPacker::transformSimilarityToLaplacian(MatrixXd& graph, MatrixXd& degrees) const
{
	unsigned N = graph.rows();
	VectorXd ones = VectorXd::Ones(N);
	VectorXd dvec = graph*ones;
	degrees = dvec.asDiagonal();

	graph *= -1;
	graph += degrees;
}

//append a single cluster
void SpectralPacker::createCluster(const DataViewer* dv, const vector<EigenInd>& vals,
		unsigned start, unsigned end, vector<Cluster>& clusters) const
{
	if(start >= end) return;
	unsigned N = end-start;

	const MappableOctTree *mivs[N];
	const MappableOctTree *msvs[N];

	vector<unsigned> indices(N);

	for(unsigned i = 0; i < N; i++)
	{
		unsigned index = vals[start+i].index;
		indices[i] = index;
		mivs[i] = dv->getMIV(index);
		msvs[i] = dv->getMSV(index);
	}

	clusters.push_back(Cluster(indices, mivs, msvs));
}

//given the ordering in vals between start and end, create packed clusters
//The goal is to have all clusters as dense as possible
void SpectralPacker::createDenseClusters(const DataViewer* dv, const vector<EigenInd>& vals, unsigned start, unsigned end, SolverPtr solver, vector<Cluster>& clusters) const
{
	if(start >= end)
		return;
	unsigned sz = end - start;

	if(sz % packSize == 0) //evenly divided, go for it
	{
		for(unsigned i = start; i < end; i += packSize)
		{
			createCluster(dv, vals, i, i+packSize, clusters);
		}
	}
	else //peel off the largest chunk we can to get balanced clusters
	{
		unsigned d = ceil(sz / (double)packSize);
		unsigned min = sz / d;

		createCluster(dv, vals, start, start+min, clusters);
		createDenseClusters(dv, vals, start+min, end, solver, clusters);
	}
}


void SpectralPacker::pack(const DataViewer* dv, const vector<unsigned>& indices,
		vector<Cluster>& clusters) const
{
	//first compute the n^2 distances
	unsigned N = indices.size();
	multi_array<double, 2> distances(extents[N][N]);
	for(unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const MappableOctTree *imiv = dv->getMIV(indices[i]);
		const MappableOctTree *imsv = dv->getMSV(indices[i]);

		for(unsigned j = 0; j < i; j++)
		{
			const MappableOctTree *jmiv = dv->getMIV(indices[j]);
			const MappableOctTree *jmsv = dv->getMSV(indices[j]);
			distances[i][j] = distances[j][i] = shapeDistance(imiv, imsv, jmiv, jmsv);
		}
	}

	//generate similarity graph
	MatrixXd graph(N,N);
	transformDistancesToSimilarity(distances, graph);
	//compute unnormalized graph laplacian
	MatrixXd D(N,N);
	transformSimilarityToLaplacian(graph, D);
	//find the second eigenvectors (first should be unit)
	SolverPtr solver;
	if(useNormalizedLaplacian)
		solver = SolverPtr(new GeneralizedSelfAdjointEigenSolver<MatrixXd>(graph, D));
	else
		solver = SolverPtr(new SelfAdjointEigenSolver<MatrixXd>(graph));

	//associate coordinates of eigenvector with indices and sort
	vector<EigenInd> vals(N);
	for(unsigned i = 0; i < N; i++)
	{
		vals[i] = EigenInd(i, indices[i], (double)solver->eigenvectors()(i,1));
	}

	//pack as densely as possible
	sort(vals.begin(), vals.end());
	createDenseClusters(dv, vals, 0, vals.size(), solver, clusters);
}
