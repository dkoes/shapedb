/*
 * Molecule.h
 *
 *  Created on: Oct 14, 2011
 *      Author: dkoes
 *
 *  This represents a molecule object.  Is supports intersection
 *  testing, file writing, and iteration over an input.
 */



#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <iostream>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>

#include "MolSphere.h"
#include "Cube.h"

using namespace OpenBabel;
using namespace std;

class MolIterator;

//interface between molecular shape and the gss tree
//defines interator type, intersection, and file create/write
class Molecule
{
	//a little icky, make mutable so intersect can bump intersecting spheres
	//to the front of the list; this ends up being more efficient than filtering
	//all the intersecting spheres
	mutable vector<MolSphere> spheres;
	OBMol mol;
public:
	typedef MolIterator iterator;

	Molecule() {}
	~Molecule() {}

	Molecule(const vector<MolSphere>& sph): spheres(sph) {}

	Molecule(const char *data)
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

	void set(const OBMol& m, float adjust=0)
	{
		mol = m;
		spheres.clear();
		spheres.reserve(mol.NumAtoms());
		for(OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
		{
			OBAtom* atom = *aitr;
			spheres.push_back(MolSphere(atom->x(), atom->y(), atom->z(), etab.GetVdwRad(atom->GetAtomicNum())+adjust));
		}
	}

	OBMol& getMol()
	{
		return mol;
	}

	bool intersects(const Cube& cube) const
	{
		for(unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			if(spheres[i].intersectsCube(cube))
			{
				swap(spheres[i],spheres[0]); //locality optimization
				return true;
			}
		}
		return false;
	}

	bool containedIn(const Cube& cube) const
	{
		for(unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			if(spheres[i].containedInCube(cube))
			{
				return true;
			}
		}
		return false;
	}

	//just return the parts of the molecule that intersect with the cube
	Molecule spliceIntersection(const Cube& cube) const
	{
		Molecule ret;
		for(unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			if(spheres[i].intersectsCube(cube))
				ret.spheres.push_back(spheres[i]);
		}
		return ret;
	}

	//return true if no shape
	bool empty() const
	{
		return spheres.size() == 0;
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
class MolIterator
{
	OBConversion inconv;
	OBMol mol;
	Molecule currmolecule;
	bool valid;
	float adjust; //chagne radii

	void readOne()
	{
		valid = inconv.Read(&mol);
		currmolecule.set(mol, adjust);
	}
public:
	MolIterator(const string& fname, float adj=0): valid(true), adjust(adj)
	{
		inconv.SetInFormat(inconv.FormatFromExt(fname));
		valid = inconv.ReadFile(&mol, fname);
		currmolecule.set(mol, adjust);
	}

	//validity check
	operator bool() const { return valid; }

	//current mol
	const Molecule& operator*() const { return currmolecule; }

	void operator++() { readOne(); }

};

#endif /* MOLECULE_H_ */
