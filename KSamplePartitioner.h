/*
 * KSamplePartitioner.h
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 */

#ifndef KSAMPLEPARTITIONER_H_
#define KSAMPLEPARTITIONER_H_

#include "GSSTreeCreator.h"
#include "MappableOctTree.h"

class KSamplePartitioner: public TopDownPartitioner
{
	unsigned kcenters; //number of centers
	unsigned ksamples; //sample ksamples*kcenters

	void kCluster(const vector<unsigned>& indices, vector< vector<unsigned> >& clusters);

	MappableOctTree* computeMIV(const vector<unsigned>& ind) const;
	MappableOctTree* computeMSV(const vector<unsigned>& ind) const;

	virtual TopDownPartitioner* create(const DataViewer* dv, vector<unsigned>& ind) const;

public:
	KSamplePartitioner(): kcenters(8), ksamples(10) {}
	KSamplePartitioner(unsigned kc, unsigned ks): kcenters(kc), ksamples(ks) {}
	KSamplePartitioner(const DataViewer *dv, unsigned kc, unsigned ks): TopDownPartitioner(dv), kcenters(kc), ksamples(ks) {}
	~KSamplePartitioner() {}

	virtual TopDownPartitioner* create(const DataViewer* dv) const;
	virtual void partition(vector<TopDownPartitioner*>& parts);
};

#endif /* KSAMPLEPARTITIONER_H_ */
