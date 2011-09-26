/*
 * main.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes

 The goal of this project is to efficiently store molecular shapes within
 a GSS-tree index for fast (and exact) similarity searching
 http://dx.doi.org/10.1145/304181.304219

 I have severe doubts that this will be able to scale the way I need it too,
 but am going to give it a shot.

 Since I intend to use this mostly as a library, the input of the primary
 routines is a list of spheres with radii already provided with a fixed
 coordinate system.  (But main will take in molecules).
 */

#define BOOST_FILESYSTEM_VERSION 2

#include <iostream>
#include <string>
#include "boost/filesystem.hpp"
#include "boost/array.hpp"
#include "openeye.h"
#include "oechem.h"
#include "oesystem.h"

#include "CommandLine2/CommandLine.h"
#include "GSSTree.h"

using namespace OEChem;
using namespace OESystem;
using namespace std;
using namespace boost;

enum CommandEnum
{
	Create, NNSearch, DCSearch, CreateSearch
};

cl::opt<CommandEnum>
		Command(
				cl::desc("Operation to perform:"),
				cl::Required,
				cl::values(clEnumVal(Create, "Create a molecule shape index."),
				clEnumVal(NNSearch, "Nearest neighbor search")						,
				clEnumVal(DCSearch, "Distance constraint search"),
				clEnumVal(CreateSearch,"In memory database creation and search for testing"),
				clEnumValEnd) );

cl::opt<string> Input("in", cl::desc("Input file"));
cl::opt<string> Output("out", cl::desc("Output file"));
cl::opt<string> Database("db", cl::desc("Database file"));

cl::opt<double> LessDist("less",
		cl::desc("Distance to reduce query mol by for constraint search (default 1A)."), cl::init(1.0));
cl::opt<double> MoreDist("more",
		cl::desc("Distance to increase query mol by for constraint search (default 1A)."), cl::init(1.0));
cl::opt<double> Resolution("resolution", cl::desc("Best resolution for shape database creation."),cl::init(1.0));

cl::opt<string> IncludeMol("incmol", cl::desc("Molecule to use for minimum included volume"));
cl::opt<string> ExcludeMol("exmol", cl::desc("Molecule to use for excluded volume"));

static void spherizeMol(OEMol& mol, vector<MolSphere>& spheres)
{
	OEAssignBondiVdWRadii(mol);
	spheres.clear();
	spheres.reserve(mol.NumAtoms());
	for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
	{
		float xyz[3];
		mol.GetCoords(atom, xyz);
		spheres.push_back(
				MolSphere(xyz[0], xyz[1], xyz[2], atom->GetRadius()));
	}

}

static void outputRes(ostream& out, const vector<MolSphere>& res)
{
	//just output in xyz w/o radii
	out << res.size();
	out << "\nTestShapeOutput\n";
	for (unsigned i = 0, n = res.size(); i < n; i++)
	{
		out << "C " << res[i].x << " " << res[i].y << " "
				<< res[i].z << "\n";
	}
}

int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);

	switch (Command)
	{
	case Create:
	case CreateSearch:
	{
		//read in all the molecules and calculate the max bounding box
		oemolistream inmols(Input);
		OEMol mol;
		float ctr[3], ext[3];
	    array<float,6> box = { {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MIN,FLT_MIN,FLT_MIN} };

	    vector<vector<MolSphere> > allmols;

	    while (OEReadMolecule(inmols, mol))
		{
	    	OEAssignBondiVdWRadii(mol);

			for (OEIter<OEConfBase> conf = mol.GetConfs(); conf; ++conf)
			{
				OEGetCenterAndExtents(conf, ctr, ext);
				box[0] = min(box[0], ctr[0] - ext[0]);
				box[1] = max(box[3], ctr[0] + ext[0]);
				box[2] = min(box[1], ctr[1] - ext[1]);
				box[3] = max(box[4], ctr[1] + ext[1]);
				box[4] = min(box[2], ctr[2] - ext[2]);
				box[5] = max(box[5], ctr[2] + ext[2]);

				vector<MolSphere> spheres;
				spheres.reserve(conf->NumAtoms());
				for(OEIter<OEAtomBase> atom = conf->GetAtoms(); atom; ++atom)
				{
					float xyz[3];
					conf->GetCoords(atom, xyz);
					spheres.push_back(MolSphere(xyz[0], xyz[1], xyz[2], atom->GetRadius()));
				}
				allmols.push_back(spheres);
			}

		}

	    //create gss tree
	    GSSTree gss(box, Resolution);
	    gss.load(allmols);

	    if(Command == Create)
	    {
	    	ofstream out(Output.c_str());
	    	gss.write(out);
	    }
	    else //create and search (in mem test)
	    {
	    	inmols.rewind();
			ofstream out(Output.c_str());
	    	while(OEReadMolecule(inmols, mol))
	    	{
	    		vector<MolSphere> spheres;
				spherizeMol(mol, spheres);
				vector<MolSphere> res;
				gss.nn_search(spheres, res);
				//just output in xyz w/o radii
				outputRes(out, res);
			}

	    }
	}
		break;
	case NNSearch:
	{
		//read in database
		GSSTree gss;
		ifstream dbfile(Database.c_str());
		if(!dbfile)
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		gss.read(dbfile);

		//read query molecule(s)
		oemolistream inmols(Input);
		OEMol mol;
		ofstream out(Output.c_str());
    	while(OEReadMolecule(inmols, mol))
    	{
    		vector<MolSphere> spheres;
			spherizeMol(mol, spheres);
			vector<MolSphere> res;
			gss.nn_search(spheres, res);
			//just output in xyz w/o radii
			outputRes(out, res);
		}

	}
	break;
	case DCSearch:
	{
		//read in database
		GSSTree gss;
		ifstream dbfile(Database.c_str());
		if(!dbfile)
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		gss.read(dbfile);
		ofstream out(Output.c_str());

		//read query molecule(s)
		if(IncludeMol.size() > 0 || ExcludeMol.size() > 0)
		{
			//use explicit volumes
			oemolistream inmol(IncludeMol);
			oemolistream exmol(ExcludeMol);
			if(!inmol && !exmol)
			{
				cerr << "Error reading inclusive/exclusive molecules\n";
				exit(-1);
			}
			vector<MolSphere> insphere, exsphere;
			if(inmol)
			{
				OEMol mol;
				OEReadMolecule(inmol, mol);
				spherizeMol(mol, insphere);
				//adjust radii
				for (unsigned i = 0, n = insphere.size(); i < n; i++)
				{
					insphere[i].incrementRadius(-LessDist);
				}
			}
			if(exmol)
			{
				OEMol mol;
				OEReadMolecule(exmol, mol);
				spherizeMol(mol, exsphere);
				//adjust radii
				for (unsigned i = 0, n = exsphere.size(); i < n; i++)
				{
					exsphere[i].incrementRadius(-MoreDist);
				}
			}

			vector<vector<MolSphere> > res;
			gss.inex_search(insphere, exsphere, res);
			//just output in xyz w/o radii
			for (unsigned i = 0, n = res.size(); i < n; i++)
				outputRes(out, res[i]);
		}
		else // range from single molecules
		{
			oemolistream inmols(Input);
			OEMol mol;
			while (OEReadMolecule(inmols, mol))
			{
				vector<MolSphere> littlespheres, bigspheres;
				spherizeMol(mol, littlespheres);
				bigspheres = littlespheres;

				//adjust radii
				for (unsigned i = 0, n = littlespheres.size(); i < n; i++)
				{
					littlespheres[i].incrementRadius(-LessDist);
					bigspheres[i].incrementRadius(MoreDist);
				}

				vector<vector<MolSphere> > res;
				gss.dc_search(littlespheres, bigspheres, res);
				//just output in xyz w/o radii
				for (unsigned i = 0, n = res.size(); i < n; i++)
					outputRes(out, res[i]);
			}
		}
	}
		break;
	}
	return 0;
}
