/*
 * MGrid.cpp
 *
 *  Created on: Nov 21, 2011
 *      Author: dkoes
 */

#include "MGrid.h"


/*
 * Molecule.cpp
 *
 *  Created on: Oct 14, 2011
 *      Author: dkoes
 */

#include "Molecule.h"

void MGrid::gridToPoint(unsigned long g, double& x, double& y, double& z) const
{
	double res = resolution;
	unsigned length = dimension/res;
	z = g % length;
	g /= length;
	y = g % length;
	g /= length;
	x = g;

	x *= res;
	y *= res;
	z *= res;

	x -= dimension/2;
	y -= dimension/2;
	z -= dimension/2;

	x += res/2.0;
	y += res/2.0;
	z += res/2.0;
}


int MGrid::pointToGrid(double x, double y, double z) const
{
	double dim = dimension;
	double res = resolution;
	unsigned len = dim/res;
	x+= dim/2;
	y+= dim/2;
	z += dim/2;

	if(x < 0) return -1;
	if(x >= dim) return -1;
	if(y < 0) return -1;
	if(y >= dim) return -1;
	if(z < 0) return -1;
	if(z >= dim) return -1;

	x /= res;
	y /= res;
	z /= res;

	unsigned ret = x;
	ret = ret*len + y;
	ret = ret*len + z;
	return ret;
}


//return the "radius" (half length) of a chord a distance d from the
//center of a circle of radius r
static double chordRadius(double r, double d)
{
	if (d > r)
		return 0;
	return sqrt(r * r - d * d);
}

//mark a line length 2*r in the grid that's a single chord of a sphere at x,y, z
void MGrid::markZChord(double x, double y, double z, double r)
{
	int start = pointToGrid(x, y, z - r);
	int end = pointToGrid(x, y, z + r);
	if(start < 0 || end < 0) return;
	grid.set_range(start, end, true);
}

//makr a circle in grid centered at x,y,z with radius r
void MGrid::markYZCircle(double x, double y, double z, double r)
{
	for (double d = 0; d <= r; d += resolution)
	{
		double cr = chordRadius(r, d);
		markZChord(x, y + d, z, cr);
		if (d != 0)
			markZChord(x, y - d, z, cr);
	}
}

void MGrid::markXYZSphere(double x, double y, double z, double r)
{
	//mark all yz circles
	for (double d = 0; d <= r; d += resolution)
	{
		double cr = chordRadius(r, d);
		markYZCircle(x + d, y, z, cr);
		if (d != 0)
			markYZCircle(x - d, y, z, cr);
	}
}

//this downsamples from from
void MGrid::copyFrom(const MGrid& from)
{
	grid.clear();
	bvect::enumerator en = from.grid.first();
	bvect::enumerator en_end = from.grid.end();

	while (en < en_end)
	{
		unsigned g = *en;
		double x,y,z;
		from.gridToPoint(g, x,y,z);
		grid.set(pointToGrid(x,y,z));
		++en;
	}
}

