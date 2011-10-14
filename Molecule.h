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
#include "openeye.h"
#include "oechem.h"
#include "oesystem.h"
#include "MolSphere.h"
#include "Cube.h"

using namespace OEChem;
using namespace OESystem;
using namespace std;

class MolIterator;

//interface between molecular shape and the gss tree
//defines interator type, intersection, and file create/write
class Molecule
{
	vector<MolSphere> spheres;
	OEMol mol;
public:
	typedef MolIterator iterator;

	Molecule() {}
	~Molecule() {}

	Molecule(const char *data)
	{
		oemolistream inm;
		inm.SetFormat(OEFormat::OEB);
		inm.Setgz(true);
		inm.openstring(data);
		OEReadMolecule(inm, mol);

		//don't bother with spheres
	}

	void set(const OEMol& m)
	{
		mol = m;
		spheres.clear();
		spheres.reserve(mol.NumAtoms());
		for(OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
		{
			float xyz[3];
			mol.GetCoords(atom, xyz);
			spheres.push_back(MolSphere(xyz[0], xyz[1], xyz[2], atom->GetRadius()));
		}
	}

	bool intersects(const Cube& cube)
	{
		for(unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			if(spheres[i].intersectsCube(cube))
				return true;
		}
		return false;
	}

	void write(ostream& out) const
	{
		oemolostream ostr;
		ostr.openstring(true);
		ostr.SetFormat(OEFormat::OEB);
		OEWriteMolecule(ostr, mol);
		const string& mstr = ostr.GetString();
		unsigned n = mstr.length();

		out.write(mstr.c_str(), n);
	}
};

// use openeye to read in from a file
class MolIterator
{
	oemolistream inmols;
	OEMol mol;
	Molecule currmolecule;
	bool valid;
	OEIter<OEConfBase> conf;

	void readOne()
	{
		if(!conf) //read next molecule
		{
			if(!OEReadMolecule(inmols, mol))
			{
				valid = false;
				return;
			}
			OEAssignBondiVdWRadii(mol);
			conf = mol.GetConfs();
		}
		//conf is now valid
		currmolecule.set(*conf);
		++conf;
	}
public:
	MolIterator(const string& fname): inmols(fname.c_str()), valid(true)
	{
		oemolistream inmols(fname);
		readOne();
	}

	//validity check
	operator bool() const { return valid; }

	//current mol
	const Molecule& operator*() const { return currmolecule; }

	void operator++() { readOne(); }

};

#endif /* MOLECULE_H_ */
