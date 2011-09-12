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

	MolSphere(float X, float Y, float Z, float R): x(X), y(Y), z(Z), r(R) {}
};



#endif /* MOLSPHERE_H_ */
