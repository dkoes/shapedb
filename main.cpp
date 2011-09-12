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

cl::opt<string> Input("in", cl::desc("Input file"), cl::Required);
cl::opt<string> Output("out", cl::desc("Output file"), cl::Required);

cl::opt<double> Distance("distance",
		cl::desc("Distance for constraint earch (default 1A)."), cl::init(1.0));
cl::opt<double> Resolution("resolution", cl::desc("Best resolution for shape database creation."),cl::init(1.0));

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

	    while (OEReadMolecule(inmols, mol))
		{
			for (OEIter<OEConfBase> conf = mol.GetConfs(); conf; ++conf)
			{
				OEGetCenterAndExtents(conf, ctr, ext);
				box[0] = min(box[0], ctr[0] - ext[0]);
				box[1] = max(box[3], ctr[0] + ext[0]);
				box[2] = min(box[1], ctr[1] - ext[1]);
				box[3] = max(box[4], ctr[1] + ext[1]);
				box[4] = min(box[2], ctr[2] - ext[2]);
				box[5] = max(box[5], ctr[2] + ext[2]);
			}

		}

	    //create gss tree
	    GSSTree gss(box, Resolution);

		//reread and generate molspheres to add to GSSTree
		inmols.rewind();
	    while (OEReadMolecule(inmols, mol))
		{
	    	OEAssignBondiVdWRadii(mol);
			for (OEIter<OEConfBase> conf = mol.GetConfs(); conf; ++conf)
			{
				vector<MolSphere> spheres;
				spheres.reserve(conf->NumAtoms());
				for(OEIter<OEAtomBase> atom = conf->GetAtoms(); atom; ++atom)
				{
					float xyz[3];
					conf->GetCoords(atom, xyz);
					spheres.push_back(MolSphere(xyz[0], xyz[1], xyz[2], atom->GetRadius()));
				}
				gss.add(spheres);
			}
		}

	    if(Command == Create)
	    {
	    	gss.write(filesystem::path(Output));
	    }
	    else //create and search (in mem test)
	    {
	    	inmols.rewind();
	    	//grab first mol
	    	OEReadMolecule(inmols, mol);
	    	vector<MolSphere> spheres;
	    	spheres.reserve(mol.NumAtoms());
	    	for(OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
	    	{
	    		float xyz[3];
	    		mol.GetCoords(atom, xyz);
	    		spheres.push_back(MolSphere(xyz[0], xyz[1], xyz[2], atom->GetRadius()));
	    	}
	    	vector<MolSphere> res;
	    	gss.nn_search(spheres, res);
	    	//just output in xyz w/o radii
	    	ofstream out(Output.c_str());
	    	out << res.size();
	    	out << "\nTestShapeOutput\n";
	    	for(unsigned i = 0, n = res.size(); i < n; i++)
	    	{
	    		out << "C " << res[i].x << " " << res[i].y << " " << res[i].z << "\n";
	    	}
	    }
	}
		break;
	case NNSearch:

	case DCSearch:
		break;
	}
	return 0;
}
