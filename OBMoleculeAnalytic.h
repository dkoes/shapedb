/*
 * OBMoleculeAnalytic.h
 *
 *  Created on: Feb 21, 2012
 *      Author: dkoes
 *
 *  This represents a molecule object.  Is supports intersection
 *  testing, file writing, and iteration over an input.
 *
 *  This version does not use grids.
 */

#ifndef OBMOLECULEA_H_
#define OBMOLECULEA_H_

#include <iostream>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <cmath>
#include "MolSphere.h"
#include "Cube.h"
#include "MGrid.h"

using namespace OpenBabel;
using namespace std;

class OBAMolIterator;
class OBAMolOutput;

/* given three close spheres and probe radius, does the necessary
 * calculations to determine if a point is inaccessible by the probe
 */
class TriChecker
{
	static float distSq(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		float X = (x2-x1);
		X *= X;
		float Y = (y2-y1);
		Y *= Y;
		float Z = (z2-z1);
		Z *= Z;

		return X+Y+Z;
	}

public:
	TriChecker(const MolSphere& x, const MolSphere& y, const MolSphere& z, float p);

	//return true if checking is necessary
	static bool isNonTrivial(const MolSphere& a, const MolSphere& b, const MolSphere& c, float probe)
	{
		float d = distSq(a.x, a.y, a.z, b.x,b.y,b.z);
		float r = a.r + b.r + probe*2;
		r *= r;

		return  d < r;
	}

};


/* given two close spheres and a probe radius, this class does the necessary
 * calculations to determine if a point is accessible by the probe
 */
class ToroidChecker
{
	MolSphere a;
	MolSphere b;
	float probe;
	float smallSq; //shortest distance to probe
	float largeSq; //longest, from line between a and b

	float AminxSq; //shortest distance from A allowed
	float AmaxxSq; //longest from A that is acceptable
	float BminxSq; //shortest distance from B allowed
	float BmaxxSq; //longest from B that is acceptable
	float probecY; //distance to probe center from line
	float probecX; //x-component of probe center from a

	float asq;
	float bsq;

	float AcheckR;

	static float distSq(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		float X = (x2-x1);
		X *= X;
		float Y = (y2-y1);
		Y *= Y;
		float Z = (z2-z1);
		Z *= Z;

		return X+Y+Z;
	}
public:
	//return true if x and y are close enough to need a checker
	static bool isNonTrivial(const MolSphere& a, const MolSphere& b, float probe)
	{
		float d = distSq(a.x, a.y, a.z, b.x,b.y,b.z);
		float r = a.r + b.r + probe*2;
		r *= r;

		return  d < r;
	}

	ToroidChecker(const MolSphere& x, const MolSphere& y, float p): a(x), b(y), probe(p)
	{
		//precalculate distances from the line between a and b that define the curved region
		//distances
		double ap = a.r+probe;
		double bp = b.r+probe;
		double ab = sqrt(distSq(a.x,a.y,a.z,b.x,b.y,b.z));
		//angles
		double A = acos((ap*ap+ab*ab-bp*bp)/(2*ap*ab));
		double B = acos((bp*bp+ab*ab-ap*ap)/(2*bp*ab));

		//bound both A and B distances since distance is unsigned
		double mina = a.r*cos(A);
		double minb = b.r*cos(B);

		AminxSq = mina*mina;
		AmaxxSq = ab-minb;
		AmaxxSq *= AmaxxSq;

		BminxSq = minb*minb;
		BmaxxSq = ab - mina;
		BmaxxSq *= BmaxxSq;

		double maxa = a.r*sin(A);
		double maxb = b.r*sin(B);
		largeSq = max(maxa,maxb);
		largeSq *= largeSq;
		//also the projection of the probe center
		probecY = (a.r+probe)*sin(A);
		probecX = (a.r+probe)*cos(A);
		smallSq = max(probecY-probe,0.0f);
		smallSq *= smallSq;

		//calculate radius from A that includes the full toroid
		float axr = ab-minb;
		float ayr = maxb;
		asq = axr*axr+ayr*ayr;

		float bxr = ab-mina;
		float byr = maxa;
		bsq = bxr*bxr+byr*byr;
	}

	float getACheckRadius() const { return sqrt(asq); }

	//return true of x,y,z is within the gap between a and b
	bool containsPoint(float x, float y, float z) const
	{
		//find the closest point on the line between a and b
		float Ax = a.x, Ay = a.y, Az = a.z;
		float Bx = b.x, By = b.y, Bz = b.z;
		float Cx = x, Cy = y, Cz = z;

		//must be within SA radius
		float Bsq = distSq(Bx,By,Bz,Cx,Cy,Cz);
		if(Bsq > bsq)
			return false;
		float Asq = distSq(Ax,Ay,Az,Cx,Cy,Cz);
		if(Asq > asq)
			return false;

		//u is the slop
		float t1 = Ax * Ax;
		float t3 = Bx * Ax;
		float t5 = Ay * Ay;
		float t7 = By * Ay;
		float t9 = Az * Az;
		float t11 = Bz * Az;
		float t13 = t1 - Ax * Cx - t3 + Bx * Cx + t5 - Ay * Cy - t7 + By * Cy + t9 - Az * Cz - t11 + Bz * Cz;
		float t15 = Bx * Bx;
		float t17 = By * By;
		float t19 = Bz * Bz;
		float u = t13 / (t1 - 2.0 * t3 + t15 + t5 - 0.2e1 * t7 + t17 + t9 - 2.0 * t11 + t19);


		//calculate actual point on line
		float Px = Ax + u*(Bx-Ax);
		float Py = Ay + u*(By-Ay);
		float Pz = Az + u*(Bz-Az);

		//distance with A
		float Ad = distSq(Ax,Ay,Az,Px,Py,Pz);
		if(Ad < AminxSq || Ad > AmaxxSq)
			return false;
		float Bd = distSq(Bx,By,Bz,Px,Py,Pz);
		if(Bd < BminxSq || Bd > BmaxxSq)
			return false;

		//distance between C and P
		float d = distSq(Px,Py,Pz,Cx,Cy,Cz);
		if(d < smallSq) //definitely fits
		{
			return true;
		}
		if(d > largeSq)
			return false;


		//falls somewhere in the curvature region
		//compute distance from probe
		float xd = probecX-Ad;
		xd *= xd;
		float yd = probecY-d;
		yd *= yd;
		if(xd+yd < probe*probe)
			return false;
		return true;
	}
};

/* This class performs intersection testing of cubes and points within the
 * context of a set of connected spheres.  It will check the sphere radius and
 * any connected toroids and triploids.
 */
class SphereChecker
{
	MolSphere atom;
	MolSphere sphere; //sphere with adjusted check radius
	vector<ToroidChecker> toroids;
public:
	SphereChecker(const MolSphere& a): atom(a), sphere(a)
	{
	}

	friend void swap(SphereChecker& a, SphereChecker& b)
	{
		using std::swap;
		// bring in swap for built-in types

		swap(a.atom,b.atom);
		swap(a.sphere,b.sphere);
		swap(a.toroids,b.toroids);
	}

	//calculate toroids and triploid if necessary
	//update radius of sphere to max check distance
	void addNeighbor(const MolSphere& b, float probe)
	{
		if(probe > 0 && ToroidChecker::isNonTrivial(atom, b,probe))
		{
			toroids.push_back(ToroidChecker(atom, b, probe));
			float checkR = toroids.back().getACheckRadius();
			if(checkR > sphere.r)
				sphere.r = checkR;
		}
	}

	//cube possibly intersects the represented volume
	bool intersectsCube(const Cube& cube) const
	{
		return sphere.intersectsCube(cube);
	}

	//point is definitely contained in represented volume
	bool containsPoint(float x, float y, float z) const
	{
		if(atom.containsPoint(x,y,z))
			return true;

		if(sphere.containsPoint(x,y,z))
		{
			for(unsigned i = 0, n = toroids.size(); i < n; i++)
			{
				if(toroids[i].containsPoint(x,y,z))
					return true;
			}
		}

		return false;
	}

};

//interface between molecular shape and the gss tree
//defines interator type, intersection, and file create/write
class OBAMolecule
{
	//a little icky, make mutable so intersect can bump intersecting spheres
	//to the front of the list; this ends up being more efficient than filtering
	//all the intersecting spheres
	mutable vector<SphereChecker> checkers;

	OBMol mol;
	float adjustAmount;
	float probe;

	OBMol& getMol()
	{
		return mol;
	}

	void set(const OBMol& m, float dimension, float resolution, float prb =
			1.4, float adj = 0)
	{
		mol = m;
		adjustAmount = adj;
		probe = prb;
		vector<MolSphere> spheres;
		spheres.reserve(mol.NumAtoms());
		for(OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
		{
			OBAtom* atom = *aitr;
			spheres.push_back(MolSphere(atom->x(), atom->y(), atom->z(), etab.GetVdwRad(atom->GetAtomicNum())-adjustAmount));
		}

		checkers.clear();
		checkers.reserve(spheres.size());
		for(unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			checkers.push_back(SphereChecker(spheres[i]));
			for(unsigned j = i+1; j < n; j++)
			{
				checkers.back().addNeighbor(spheres[j],probe);
			}
		}
	}


	friend class OBAMolOutput; //needs to call getMol
	friend class OBAMolIterator; //calls set

public:
	typedef OBAMolIterator iterator;
	typedef OBAMolOutput molostream;

	OBAMolecule(): adjustAmount(0), probe(0) {}
	~OBAMolecule()
	{
	}


	OBAMolecule(const char *data): adjustAmount(0), probe(0)
	{
		OBConversion conv;
		conv.SetInAndOutFormats("SDF","SDF");
		conv.SetOptions("z",OBConversion::INOPTIONS);
		conv.SetOptions("z",OBConversion::GENOPTIONS);

		unsigned n = 0;
		memcpy(&n, data, sizeof(unsigned));
		const char* mdata = (const char*)data+sizeof(unsigned);

		string mstr(mdata, n);

		conv.ReadString((OBBase*)&mol, mstr);

		//don't bother with spheres
	}


	//mutate small to be smaller/larger
	void adjust(float dimension, float resolution, float probe, float adj)
	{
		set(mol, dimension, resolution, probe, adj);
	}

	bool intersects(const Cube& cube) const
	{
		for(unsigned i = 0, n = checkers.size(); i < n; i++)
		{
			if(checkers[i].intersectsCube(cube))
			{
				if(i > 0) swap(checkers[i],checkers[0]); //locality optimization

				return true;
			}
		}
		return false;
	}

	//return true if point is within object
	bool containsPoint(float x, float y, float z) const
	{
		for(unsigned i = 0, n = checkers.size(); i < n; i++)
		{
			if(checkers[i].containsPoint(x,y,z))
			{
				if(i > 0) swap(checkers[i],checkers[0]); //locality optimization
				return true;
			}
		}

		return false;
	}


	void write(ostream& out) const
	{
		OBConversion conv;
		conv.SetInAndOutFormats("SDF","SDF");
		conv.SetOptions("z",OBConversion::OUTOPTIONS);
		conv.SetOptions("z",OBConversion::GENOPTIONS);

		OBMol copy = mol; //because OB doesn't have a const version
		string mstr = conv.WriteString(&copy);

		unsigned n = mstr.length();
		out.write((char*)&n, sizeof(unsigned));
		out.write(mstr.c_str(), n);
	}


};


// use openeye to read in from a file
class OBAMolIterator
{
	OBConversion inconv;
	OBMol mol;
	OBAMolecule currmolecule;
	bool valid;
	float dimension;
	float resolution;
	float probe;
	float adjust; //chagne radii

	void readOne()
	{
		valid = inconv.Read(&mol);
		currmolecule.set(mol, dimension, resolution, probe, adjust);
	}
public:
	OBAMolIterator(const string& fname, float dim, float res, float prb = 1.4,
			float adj = 0) :
			valid(true), dimension(dim), resolution(res), probe(prb), adjust(
					adj)
	{
		inconv.SetInFormat(inconv.FormatFromExt(fname));
		valid = inconv.ReadFile(&mol, fname);
		currmolecule.set(mol, dimension, resolution, probe, adjust);
	}

	//validity check
	operator bool() const
	{
		return valid;
	}

	//current mol
	const OBAMolecule& operator*() const
	{
		return currmolecule;
	}

	void operator++()
	{
		readOne();
	}

};

//output molecules to a file
class OBAMolOutput
{
	OBConversion outconv;
	ofstream out;
public:
	//open file for output
	OBAMolOutput(const string& fname)
	{
		outconv.SetOutFormat(outconv.FormatFromExt(fname.c_str()));
		out.open(fname.c_str());
		outconv.SetOutStream(&out);
	}

	bool write(OBAMolecule& mol)
	{
		return outconv.Write(&mol.getMol());
	}
};

#endif /* OBMOLECULEA_H_ */
