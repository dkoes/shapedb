/*
 * MGrid.h
 *
 *  Created on: Nov 21, 2011
 *  A class for representing a grid of molecular data.
 *      Author: dkoes
 */

#ifndef MGRID_H_
#define MGRID_H_

#include <bm/bm.h>

typedef bm::bvector<bm::standard_allocator> bvect;

class MGrid
{

	void gridToPoint(unsigned long g, double& x, double& y, double& z) const;
	int pointToGrid(double x, double y, double z) const;
	void markZChord(double x, double y, double z, double r);
	void markYZCircle(double x, double y, double z, double r);


	double dimension;
	double resolution;
	bvect grid;
public:
	MGrid() {}
	MGrid(float d, float r): dimension(d), resolution(r) {}
	~MGrid() {}

	void markXYZSphere(double x, double y, double z, double r);
	double getResolution() const { return resolution; }
	double getDimension() const { return dimension; }

	void copyFrom(const MGrid& from); //keeps current grid dim/res, samples from from

	bool test(float x, float y, float z) const
	{
		return grid.test(pointToGrid(x,y,z));
	}
};

#endif /* MGRID_H_ */
