/*
 * Packer.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef PACKER_H_
#define PACKER_H_

#include "DataViewers.h"
//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
protected:
	unsigned packSize;
public:

	Packer(unsigned ps): packSize(ps) {}
	virtual ~Packer() {}
	virtual void pack(const DataViewer* dv, const vector<unsigned>& indices, vector<Cluster>& clusters) const = 0;
};



#endif /* PACKER_H_ */
