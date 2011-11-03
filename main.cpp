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


#include <iostream>
#include <string>
#include "GSSTypes.h"
#include "boost/array.hpp"
#include "openeye.h"
#include "oechem.h"
#include "oesystem.h"

#include "CommandLine2/CommandLine.h"
#include "GSSTreeCreator.h"
#include "GSSTreeSearcher.h"
#include "Molecule.h"
#include "KSamplePartitioner.h"
#include "FullMergePacker.h"
#include "MatcherPacker.h"
#include "SpectralPacker.h"
#include "GreedyPacker.h"
#include <boost/shared_ptr.hpp>

using namespace OEChem;
using namespace OESystem;
using namespace std;
using namespace boost;

typedef shared_ptr<Packer> PackerPtr;

enum CommandEnum
{
	Create, NNSearch, DCSearch
};

cl::opt<CommandEnum>
		Command(
				cl::desc("Operation to perform:"),
				cl::Required,
				cl::values(clEnumVal(Create, "Create a molecule shape index."),
				clEnumVal(NNSearch, "Nearest neighbor search")						,
				clEnumVal(DCSearch, "Distance constraint search"),
				clEnumValEnd) );

enum PackerEnum
{
	FullMerge, Spectral,GreedyMerge,MatchPack
};
cl::opt<PackerEnum>
		Packer(
				cl::desc("Packing algorithm:"),
				cl::values(clEnumValN(FullMerge,"full-merge", "Greedy full merge."),
						clEnumValN(GreedyMerge,"greedy-merge", "Greedy iterative merge."),
						clEnumValN(MatchPack,"match-merge", "Optimal matching merging."),
				clEnumValN(Spectral, "spectral", "Spectral packing"),
				clEnumValEnd), cl::init(MatchPack) );

cl::opt<SpectralPacker::SpectralAlgEnum>
	SpectralAlg(cl::desc("Spectral packing sub-algorithm:"),
			cl::values(clEnumValN(SpectralPacker::SortDense, "sort-dense", "Simple sort followed by dense packing"),
					clEnumValN(SpectralPacker::SortPartition, "sort-partition", "Simple sort followed by largest separator partition packing"),
					clEnumValN(SpectralPacker::ClusterFullEigen, "cluster-eigen", "Cluster eigen values using greedy packer"),
					clEnumValN(SpectralPacker::RelaxationPacking, "relax", "Cluster form relaxation values"),
					clEnumValEnd), cl::init(SpectralPacker::SortDense));

cl::opt<bool> ScanCheck("scancheck",
		cl::desc("Perform a full scan to check results"), cl::Hidden);
cl::opt<bool> ScanOnly("scanonly",
		cl::desc("Search using only a scan"), cl::Hidden);

cl::opt<string> Input("in", cl::desc("Input file"));
cl::opt<string> Output("out", cl::desc("Output file"));
cl::opt<string> Database("db", cl::desc("Database file"));

cl::opt<double> LessDist("less",
		cl::desc("Distance to reduce query mol by for constraint search (default 1A)."), cl::init(1.0));
cl::opt<double> MoreDist("more",
		cl::desc("Distance to increase query mol by for constraint search (default 1A)."), cl::init(1.0));
cl::opt<double> MaxDimension("max-dim", cl::desc("Maximum dimension."),cl::init(64));
cl::opt<double> Resolution("resolution", cl::desc("Best resolution for shape database creation."),cl::init(.5));

cl::opt<string> IncludeMol("incmol", cl::desc("Molecule to use for minimum included volume"));
cl::opt<string> ExcludeMol("exmol", cl::desc("Molecule to use for excluded volume"));

cl::opt<unsigned> KCenters("kcenters", cl::desc("number of centers for ksample-split"), cl::init(8));
cl::opt<unsigned> KSampleMult("ksamplex", cl::desc("multiplictive factor for ksampling"),cl::init(10));

cl::opt<unsigned> LeafPack("leaf-pack",
		cl::desc("Cutoff to trigger leaf packing"), cl::init(256));
cl::opt<unsigned> NodePack("node-pack",
		cl::desc("Cutoff to trigger leaf packing"), cl::init(256));

cl::opt<unsigned> Pack("pack", cl::desc("Maximum quantities per a node"), cl::init(8));

cl::opt<bool> Verbose("v", cl::desc("Verbose output"));
cl::opt<bool> UseUnnorm("use-unnorm", cl::desc("Use unnormalized laplacian in spectral packing"), cl::init(false));

cl::opt<Packer::ClusterDistance>		ClusterDist(
				cl::desc("Metric for cluster packing distance:"),
				cl::values(clEnumValN(Packer::AverageLink, "ave-dist", "Use 'average' metric between MIV/MSV representations of clusters"),
						clEnumValN(Packer::CompleteLink, "complete-dist", "Use complete linkage value between cluster members"),
						clEnumValN(Packer::SingleLink, "single-dist", "Use single linkage value between cluster members"),
				clEnumValEnd),cl::init(Packer::AverageLink) );


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


int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);

	switch (Command)
	{
	case Create:
	{
		//read in all the molecules and calculate the max bounding box
		KSamplePartitioner topdown(KCenters, KSampleMult);

		PackerPtr packer;
		switch(Packer)
		{
		case FullMerge:
			packer = PackerPtr(new FullMergePacker(Pack, ClusterDist));
			break;
		case MatchPack:
			packer = PackerPtr(new MatcherPacker(Pack, ClusterDist));
			break;
		case GreedyMerge:
			packer = PackerPtr(new GreedyPacker(Pack, ClusterDist));
			break;
		case Spectral:
			packer = PackerPtr(new SpectralPacker(Pack, SpectralAlg, !UseUnnorm));
			break;
		}


		GSSLevelCreator leveler(&topdown, packer.get(), NodePack, LeafPack);

		GSSTreeCreator creator(&leveler);

		filesystem::path dbpath(Database.c_str());
		Molecule::iterator molitr(Input);
		if(!creator.create(dbpath, molitr, MaxDimension, Resolution))
		{
			cerr << "Error creating database\n";
			exit(1);
		}

		if(Verbose)
			creator.printStats(cout);

	}
		break;
	case NNSearch:
	{
		//read in database
		ifstream dbfile(Database.c_str());
		if(!dbfile)
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		//TODO: load db

		//read query molecule(s)
		for(Molecule::iterator molitr(Input); molitr; ++molitr)
		{

		}


	}
	break;
	case DCSearch:
	{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if(!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}


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

			vector<Molecule> res;

			//search
			if(!ScanOnly)
				gss.dc_search(Molecule(insphere), Molecule(exsphere), true, res);

			if(ScanCheck || ScanOnly)
			{
				vector<Molecule> res2;
				gss.dc_scan_search(Molecule(insphere), Molecule(exsphere), true, res2);
				if(res2.size() != res.size())
				{
					cerr << "Scanning found different number\n";
				}
			}

			oemolostream out(Output.c_str());
			for (unsigned i = 0, n = res.size(); i < n; i++)
				res[i].writeMol(out);
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

				vector<Molecule > res;
				//search
				if(!ScanOnly)
					gss.dc_search(Molecule(littlespheres), Molecule(bigspheres), false, res);

				if(ScanCheck || ScanOnly)
				{
					vector<Molecule> res2;
					gss.dc_scan_search(Molecule(littlespheres), Molecule(bigspheres), false, res2);
					if(res2.size() != res.size())
					{
						cerr << "Scanning found different number\n";
					}
				}

				oemolostream out(Output.c_str());
				for (unsigned i = 0, n = res.size(); i < n; i++)
					res[i].writeMol(out);
			}
		}
	}
		break;
	}
	return 0;
}
