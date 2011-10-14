/*
 * MolSphere.h
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#ifndef MOLSPHERE_H_
#define MOLSPHERE_H_

#include "Cube.h"
//a sphere
struct MolSphere
{
	float x;
	float y;
	float z;
	float r;

	MolSphere(): x(0), y(0), z(0), r(0) {}
	MolSphere(float X, float Y, float Z, float R): x(X), y(Y), z(Z), r(R) {}

	bool operator==(const MolSphere& rhs) const {
		return x == rhs.x && y == rhs.y && z == rhs.z && r == rhs.r;
	}

	void incrementRadius(float diff)
	{
		r += diff;
		if(r < 0) r = 0;
	}

	inline float squared(float v) const
	{
		return v * v;
	}

	//thank you stack overflow for not making me think..
	bool intersectsCube(const Cube& cube) const
	{
		float dist = r * r;
		float x2 = cube.x + cube.dim;
		float y2 = cube.y + cube.dim;
		float z2 = cube.z + cube.dim;

		if (x < cube.x)
			dist -= squared(x - cube.x);
		else if (x > x2)
			dist -= squared(x - x2);

		if (y < cube.y)
			dist -= squared(y - cube.y);
		else if (y > y2)
			dist -= squared(y - y2);

		if (z < cube.z)
			dist -= squared(z - cube.z);
		else if (z > z2)
			dist -= squared(z - z2);

		return dist > 0;
	}

};



#endif /* MOLSPHERE_H_ */
