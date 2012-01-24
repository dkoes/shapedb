/*
 * OEMolecule.h
 *
 *  Created on: Jan 24, 2011
 *      Author: dkoes
 *
 *  This represents a molecule object.  Is supports intersection
 *  testing, file writing, and iteration over an input.
 *
 *  Uses openeye tools on the backend
 */

#ifndef OEMOLECULE_H_
#define OEMOLECULE_H_

#include <iostream>
#include <string>


#include <openeye.h>
#include <oechem.h>
#include <oesystem.h>
#include <oespicoli.h>

#include "MolSphere.h"
#include "Cube.h"
#include "MGrid.h"

using namespace std;
using namespace OEChem;
using namespace OESystem;
using namespace OESpicoli;

class OEMolIterator;
class OEMolOutput;

//interface between molecular shape and the gss tree
//defines interator type, intersection, and file create/write
class OEMolecule
{
	OEMol mol;
	vector<OEScalarGrid> grids;

	//downsample to a lower resolution, generating an MSV
	void downsample(const OEScalarGrid& source, OEScalarGrid& target, float dim, float newres)
	{
		float minmax[6] = {-dim,-dim,-dim,dim,dim,dim};
		target = OEScalarGrid(minmax,newres);

		for(unsigned i = 0, n = source.GetSize(); i < n; i++)
		{
			if(source[i])
			{
				float x=0, y=0, z=0;
				source.ElementToSpatialCoord(i, x, y, z);

				target(x,y,z) =1 ;
			}
		}

	}
	//create a set of grids of different resolution from spheres
	void createGrids(double dimension, double resolution, double probeRadius, double adjust)
	{
		grids.clear();
		OESurface surf;

		if(probeRadius == 0) //molecular
			OEMakeMolecularSurface(surf, mol, resolution);
		else //solvent accessible
			OEMakeAccessibleSurface(surf, mol, resolution, probeRadius);

		OEScalarGrid mgrid;
		OEMakeGridFromSurface(mgrid, surf,  OEVoxelizeMethod::Distance);

		for(unsigned i = 0, n = mgrid.GetSize(); i < n; i++)
		{
			//discretize to 0/1
			if(mgrid[i] <= -adjust)
			{
				mgrid[i] = 1;
			}
			else
				mgrid[i] = 0;
		}

		OEScalarGrid grid;
		//keep resolution, be shift to canonical grid
		downsample(mgrid, grid, dimension, resolution);


		grids.push_back(grid);
		double res = resolution * 2;
		while (res <= dimension)
		{
			downsample(grids.back(), grid, dimension, res);
			grids.push_back(grid);
			res *= 2;
		}
	}

	//return mol, called from oemoloutput
	OEMol& getMol()
	{
		return mol;
	}

	//set mol and initialize grids, called form moliterator
	void set(const OEMol& m, float dimension, float resolution, float probe =
			1.4, float adjust = 0)
	{
		mol = m;
		createGrids(dimension, resolution, probe, adjust);
	}

public:

	friend class OEMolIterator;
	friend class OEMolOutput;
	typedef OEMolIterator iterator;
	typedef OEMolOutput molostream;

	OEMolecule()
	{
	}

	~OEMolecule()
	{
	}

	void adjust(float dimension, float resolution, float probe, float adjust)
	{
		createGrids(dimension, resolution, probe, adjust);
	}

	//read in molecular data from data stream that was written with write, must call adjust to initialize grids
	OEMolecule(const char *data)
	{
		unsigned n = 0;
		memcpy(&n, data, sizeof(unsigned));
		const unsigned char* mdata = (const unsigned char*) data + sizeof(unsigned);

		oemolistream ifs;
		ifs.SetFormat(OEFormat::OEB);
		ifs.Setgz(true);
		ifs.openstring(mdata, n);

		OEReadMolecule(ifs, mol);
	}

	//write molecular data
	void write(ostream& out) const
	{
		oemolostream ofs;
		ofs.openstring();
		ofs.SetFormat(OEFormat::OEB);
		ofs.Setgz(true);

		OEWriteMolecule(ofs, mol);
		string mstr = ofs.GetString();
		unsigned n = mstr.length();
		out.write((char*) &n, sizeof(unsigned));
		out.write(mstr.c_str(), n);

	}

	bool intersects(const Cube& cube) const
	{
		/*
		bool fullcheck = false;
		float cx=0,cy=0,cz=0;
		cube.getBottomCorner(cx,cy,cz);
		float dim = cube.getDimension();
		float res = grids[0].GetSpacing();
		for(float x = cx+res/2, stopx = cx+dim; x < stopx; x+= res)
		{
			for(float y=cy+res/2, stopy = cy+dim; y < stopy; y+= res)
			{
				for(float z=cz+res/2, stopz = cz+dim; z < stopz; z+= res)
				{
					if(grids[0].IsInGrid(x,y,z) && grids[0](x,y,z))
					{
						fullcheck = true;
						break;
					}
				}
				if(fullcheck) break;
			}
			if(fullcheck) break;
		}
		return fullcheck;
*/
		bool hcheck = false;
		for (unsigned i = 0, n = grids.size(); i < n; i++)
		{
			if (grids[i].GetSpacing() == cube.getDimension())
			{
				float x, y, z;
				cube.getCenter(x, y, z);
				if(grids[i].IsInGrid(x,y,z) && grids[i](x,y,z))
					hcheck = true;
				else
					hcheck = false;
				break;
			}
		}

		return hcheck;
	}

	//return true if no shape
	bool empty() const
	{
		return mol.NumAtoms() == 0;
	}

};

// use openeye to read in from a file
class OEMolIterator
{
	oemolistream ifs;
	OEMolecule currmolecule;
	bool valid;
	float dimension;
	float resolution;
	float probe;
	float adjust; //chagne radii

	void readOne()
	{
		OEMol mol;
		valid = OEReadMolecule(ifs, mol);
		if(valid)
			currmolecule.set(mol, dimension, resolution, probe, adjust);
	}
public:
	OEMolIterator(const string& fname, float dim, float res, float prb = 1.4,
			float adj = 0) :
			valid(true), dimension(dim), resolution(res), probe(prb), adjust(
					adj)
	{
		ifs.open(fname);
		readOne();
	}

	//validity check
	operator bool() const
	{
		return valid;
	}

	//current mol
	const OEMolecule& operator*() const
	{
		return currmolecule;
	}

	void operator++()
	{
		readOne();
	}

};

class OEMolOutput
{
	oemolostream ofs;
public:
	//open file for output
	OEMolOutput(const string& fname)
	{
		ofs.open(fname);
	}

	bool write(OEMolecule& mol)
	{
		return OEWriteMolecule(ofs, mol.getMol());
	}
};

#endif /* MOLECULE_H_ */
