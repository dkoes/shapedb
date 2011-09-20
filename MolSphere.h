/*
 * MolSphere.h
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#ifndef MOLSPHERE_H_
#define MOLSPHERE_H_

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
};



#endif /* MOLSPHERE_H_ */
