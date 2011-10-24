/*
 * SpectralPacker.h
 *
 *  Created on: Oct 24, 2011
 *      Author: dkoes
 *
 *  Packs bottom up using the spectral ordering - the second eigenvector
 *  of a graph Laplacian. Enhances partitioning using additional eigenvectors.
 */

#ifndef SPECTRALPACKER_H_
#define SPECTRALPACKER_H_

#include "Packer.h"
#include "boost/multi_array.hpp"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace boost;
using namespace Eigen;


class SpectralPacker: public Packer
{
	//small container for storing and sorting values from eigen vector and index
	struct EigenInd
	{
		double value;
		unsigned index; //value from indices
		unsigned index_of_index; //position within indices, for referencing whole eigenvalue

		EigenInd(): value(0), index(0), index_of_index(0) {}
		EigenInd(unsigned i, unsigned i2, double v): value(v), index(i2), index_of_index(i) {}

		bool operator<(const EigenInd& rhs) const
		{
			return value < rhs.value;
		}
	};
	typedef shared_ptr< SelfAdjointEigenSolver<MatrixXd> > SolverPtr;

	bool useNormalizedLaplacian;
	void transformDistancesToSimilarity(const multi_array<double, 2>& distances, MatrixXd& graph) const;
	void transformSimilarityToLaplacian(MatrixXd& graph, MatrixXd& degrees) const;

	void createCluster(const DataViewer* dv, const vector<EigenInd>& vals,
			unsigned start, unsigned end, vector<Cluster>& clusters) const;

	void createDenseClusters(const DataViewer* dv, const vector<EigenInd>& vals,
			unsigned start, unsigned end, SolverPtr solver, vector<Cluster>& clusters) const;
public:
	SpectralPacker(unsigned pk, bool norm=true): Packer(pk, NotApplicable), useNormalizedLaplacian(norm) {}
	~SpectralPacker() {}

	void pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const;

};

#endif /* SPECTRALPACKER_H_ */
