/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
/*
 * MGrid.h
 *
 *  Created on: Nov 21, 2011
 *  A class for representing a grid of molecular data.
 *      Author: dkoes
 */

#ifndef MGRID_H_
#define MGRID_H_

#include <cassert>
#include <bm/bm.h>
#include <vector>
#include "Cube.h"
#include <eigen3/Eigen/Core>

typedef bm::bvector<bm::standard_allocator> bvect;

class MGrid
{

	void gridToPoint(unsigned long g, double& x, double& y, double& z) const;
	int pointToGrid(double x, double y, double z) const;
	void markZChord(double x, double y, double z, double r);
	void markYZCircle(double x, double y, double z, double r);
	void shrinkByOne();
	void growByOne();

	 void makeFace(Eigen::Vector3f pt, int which, double r,
			 vector<Eigen::Vector3f>& vertices, vector<Eigen::Vector3f>& normals, vector<int>& faces);


	double dimension; //max size, in distance unites
	double resolution; //size of each grid cube
	bvect grid;
public:

	struct Point
	{
		float x, y, z;
		Point(): x(0), y(0), z(0) {}
		Point(float X, float Y, float Z): x(X), y(Y), z(Z) {}
	};

	MGrid(): dimension(32), resolution(0.5) {}
	MGrid(float d, float r): dimension(d), resolution(r) {}
	~MGrid() {}

	void markXYZSphere(double x, double y, double z, double r);
	double getResolution() const { return resolution; }
	double getDimension() const { return dimension; }

	void copyFrom(const MGrid& from); //keeps current grid dim/res, samples from from

	bool test(float x, float y, float z) const
	{
		int pt = pointToGrid(x,y,z);
		if(pt >= 0)
			return grid.test(pt);
		else
			return true; //for exposed point testing and shrinking this makes the most sense
	}

	bool inGrid(float x, float y, float z) const
	{
		return pointToGrid(x,y,z) >= 0;
	}

	//test for x,y,z - match signature expected from oct tree creation
	bool containsPoint(float x, float y, float z) const
	{
		return test(x,y,z);
	}

	//return true if possibly intersects cube - by always returning true
	//explores every grid point with contains point resulting in some
	//unnecessary oct-tree creation
	bool intersects(const Cube& cube) const
	{
		return true;
	}

	void setPoint(float x, float y, float z);
	void makeSurface(const MGrid& sagrid, double probe);

	bool isInteriorPoint(float x, float y, float z) const;
	bool isExposedPoint(float x, float y, float z) const;
	bool isSolitaryPoint(float x, float y, float z) const;

	//return true if at least part of the sphere fits in the grid
	bool sphereInGrid(float x, float y, float z, float r) const;
	void shrink(double amount);
	void grow(double amount);

	void operator&=(MGrid& rhs)
	{
		assert(dimension == rhs.dimension && resolution == rhs.resolution);
		grid &= rhs.grid;
	}

	void operator|=(MGrid& rhs)
	{
		assert(dimension == rhs.dimension && resolution == rhs.resolution);
		grid |= rhs.grid;
	}

	void clear()
	{
		grid.clear();
	}

	unsigned numSet() const
	{
		return grid.count();
	}

	void getSetPoints(std::vector<Point>& points) const;

	void makeMesh(vector<Eigen::Vector3f>& vertices, vector<Eigen::Vector3f>& normals, vector<int>& faces);
};

#endif /* MGRID_H_ */
