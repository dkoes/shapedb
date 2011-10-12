/*
 * OctTree.h
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 *
 *      An octtree representation of a molecular shape.
 *      Requires a predefined bounding box and resolution.
 *      Voxelizes a molecule represented as a collection of spheres.
 *
 *      This is an abstract class and a factory class so we can
 *      evaluate different data structures.
 */

#ifndef OCTTREE_H_
#define OCTTREE_H_

#include "MolSphere.h"
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include "Cube.h"

using namespace std;


class OctTree
{

public:
	OctTree()
	{
	}

	virtual ~OctTree();

	virtual OctTree* clone() const = 0;
	//invert filled and unfilled
	virtual void invert() = 0;

	//mogrifying intersection
	virtual bool intersect(const OctTree* rhs) = 0;
	//mogrifying union
	virtual bool unionWith(const OctTree* rhs) = 0;

	//volume calculations that don't require creating a tmp tree
	virtual float intersectVolume(const OctTree * rhs) const = 0;
	virtual float unionVolume(const OctTree *rhs) const = 0;
	virtual float subtractVolume(const OctTree *rhs) const
	{
		OctTree *tmp = clone();
		tmp->invert();
		tmp->intersect(rhs);
		float ret = tmp->volume();
		delete tmp;
		return ret;
	}

	//return total volume contained in octtree
	virtual float volume() const = 0;

	//return number of leaves
	virtual unsigned leaves() const = 0;

	virtual void clear() = 0;
	virtual void fill() = 0;

	virtual void write(ostream& out) const = 0;
	virtual void read(istream& in) = 0;

	virtual unsigned getOctantPattern(const vector<unsigned>& coord, bool MSV) const = 0;

	virtual float relativeVolumeDistance(const OctTree * B) const
	{
		return 1 - intersectVolume(B)/unionVolume(B);
	}

	virtual float absoluteVolumeDistance(const OctTree * B) const
	{
		return unionVolume(B)-intersectVolume(B);
	}

	virtual float hausdorffDistance(const OctTree* B) const = 0;

	virtual bool containedIn(const OctTree *larger) const
	{
		return intersectVolume(larger) == volume();
	}

	virtual float percentOverlapVolume(const OctTree *thisMSV, const OctTree *rightMIV, const OctTree *rightMSV) const
	{
		OctTree *tmp1 = clone();
		tmp1->unionWith(rightMIV);
		tmp1->invert();

		OctTree *tmp2 = thisMSV->clone();
		tmp2->intersect(rightMSV);

		tmp1->intersect(tmp2); //overlap

		OctTree *tmp3 = clone();
		tmp3->intersect(rightMIV);
		tmp3->invert();

		OctTree *tmp4 = thisMSV->clone();
		tmp4->unionWith(rightMSV);
		tmp4->intersect(tmp3); //anylap

		float ret;
		if(tmp4->volume() == 0)
			ret = 0;
		else
			ret = 1-tmp1->volume()/tmp4->volume();
		delete tmp1;
		delete tmp2;
		delete tmp3;
		delete tmp4;

		return ret;

	}


	//see how much the volume inbetween MIV and MSV will increase
	virtual float inbetweenVolume(const OctTree *thisMSV, const OctTree *rightMIV, const OctTree *rightMSV) const
	{
		OctTree *tmp1 = clone();
		tmp1->intersect(rightMIV); //new miv

		OctTree *tmp2 = thisMSV->clone();
		tmp2->unionWith(rightMSV); //new msv

		float newvol = tmp1->subtractVolume(tmp2);
		float old1 = subtractVolume(thisMSV);
		float old2 = rightMIV->subtractVolume(rightMSV);

		float ret = 1-(old1+old2)/(2*newvol);
		ret = (newvol-old1) + (newvol-old2);
		delete tmp1;
		delete tmp2;

		return ret;
	}
};



#endif /* OCTTREE_H_ */
