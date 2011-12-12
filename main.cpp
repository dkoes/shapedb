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

#include "CommandLine2/CommandLine.h"
#include "GSSTreeCreator.h"
#include "GSSTreeSearcher.h"
#include "Molecule.h"
#include "KSamplePartitioner.h"
#include "packers/Packers.h"
#include <boost/shared_ptr.hpp>
#include "Timer.h"

using namespace std;
using namespace boost;

typedef shared_ptr<Packer> PackerPtr;

enum CommandEnum
{
	Create, NNSearch, DCSearch, MolGrid, BatchSearch, BatchDB
};

cl::opt<CommandEnum> Command(
		cl::desc("Operation to perform:"),
		cl::Required,
		cl::values(clEnumVal(Create, "Create a molecule shape index."),
				clEnumVal(NNSearch, "Nearest neighbor search"),
				clEnumVal(DCSearch, "Distance constraint search"),
				clEnumVal(MolGrid, "Generate molecule grid and debug output"),
				clEnumVal(BatchSearch, "Read in a jobs file for batch processing"),
				clEnumVal(BatchDB, "Read in a jobs file for batch processing of a list of directories"),
				clEnumValEnd));

enum PackerEnum
{
	FullMerge, Spectral, GreedyMerge, MatchPack
};
cl::opt<PackerEnum> PackerChoice(
		cl::desc("Packing algorithm:"),
		cl::values(
				clEnumValN(FullMerge,"full-merge", "Greedy full merge."),
				clEnumValN(GreedyMerge,"greedy-merge", "Greedy iterative merge."),
				clEnumValN(MatchPack,"match-merge", "Optimal matching merging."),
				clEnumValN(Spectral, "spectral", "Spectral packing"),
				clEnumValEnd), cl::init(MatchPack));

cl::opt<bool> QuadPack("quad-pack", cl::desc("Use quad packing"), cl::init(false));
cl::opt<unsigned> K("k", cl::desc("k nearest neighbors to find for NNSearch"),
		cl::init(1));
cl::opt<unsigned> Knn("knn", cl::desc("K for knn graph creation"), cl::init(8));
cl::opt<unsigned> Sentinals("sentinals",
		cl::desc("Number of sentinals for knn initialization (zero random)"),
		cl::init(32));
cl::opt<SpectralPacker::SpectralAlgEnum> SpectralAlg(
		cl::desc("Spectral packing sub-algorithm:"),
		cl::values(
				clEnumValN(SpectralPacker::SortDense, "sort-dense", "Simple sort followed by dense packing"),
				clEnumValN(SpectralPacker::SortPartition, "sort-partition", "Simple sort followed by largest separator partition packing"),
				clEnumValN(SpectralPacker::ClusterFullEigen, "cluster-eigen", "Cluster eigen values using greedy packer"),
				clEnumValN(SpectralPacker::RelaxationPacking, "relax", "Cluster form relaxation values"),
				clEnumValEnd), cl::init(SpectralPacker::SortDense));

cl::opt<bool> ScanCheck("scancheck",
		cl::desc("Perform a full scan to check results"), cl::Hidden);
cl::opt<bool> ScanOnly("scanonly", cl::desc("Search using only a scan"),
		cl::Hidden);

cl::opt<string> Input("in", cl::desc("Input file"));
cl::opt<string> Output("out", cl::desc("Output file"));
cl::opt<string> Database("db", cl::desc("Database file"));

cl::opt<double> LessDist(
		"less",
		cl::desc(
				"Distance to reduce query mol by for constraint search (default 1A)."),
		cl::init(1.0));
cl::opt<double> MoreDist(
		"more",
		cl::desc(
				"Distance to increase query mol by for constraint search (default 1A)."),
		cl::init(1.0));

cl::opt<double> ProbeRadius("probe-radius",
		cl::desc("Radius of water probe for SA calculations"), cl::init(0));

cl::opt<double> MaxDimension("max-dim", cl::desc("Maximum dimension."),
		cl::init(64));
cl::opt<double> Resolution("resolution",
		cl::desc("Best resolution for shape database creation."), cl::init(.5));

cl::opt<string> IncludeMol("incmol",
		cl::desc("Molecule to use for minimum included volume"));
cl::opt<string> ExcludeMol("exmol",
		cl::desc("Molecule to use for excluded volume"));

cl::opt<unsigned> KCenters("kcenters",
		cl::desc("number of centers for ksample-split"), cl::init(8));
cl::opt<unsigned> KSampleMult("ksamplex",
		cl::desc("multiplictive factor for ksampling"), cl::init(10));

cl::opt<unsigned> LeafPack("leaf-pack",
		cl::desc("Cutoff to trigger leaf packing"), cl::init(100000));
cl::opt<unsigned> NodePack("node-pack",
		cl::desc("Cutoff to trigger leaf packing"), cl::init(100000));

cl::opt<unsigned> Pack("pack", cl::desc("Maximum quantities per a node"),
		cl::init(8));

cl::opt<bool> Verbose("v", cl::desc("Verbose output"));
cl::opt<bool> UseUnnorm("use-unnorm",
		cl::desc("Use unnormalized laplacian in spectral packing"),
		cl::init(false));

cl::opt<Packer::ClusterDistance> ClusterDist(
		cl::desc("Metric for cluster packing distance:"),
		cl::values(
				clEnumValN(Packer::AverageLink, "ave-dist", "Use 'average' metric between MIV/MSV representations of clusters"),
				clEnumValN(Packer::CompleteLink, "complete-dist", "Use complete linkage value between cluster members"),
				clEnumValN(Packer::SingleLink, "single-dist", "Use single linkage value between cluster members"),
				clEnumValN(Packer::TotalLink, "total-dist", "Use total (sum) linkage value between cluster members"),
				clEnumValEnd), cl::init(Packer::AverageLink));

cl::opt<DistanceFunction> ShapeDist(
		cl::desc("Metric for distance between shapes:"),
		cl::values(
				clEnumValN(RelativeVolume,"rel-volume","Relative volume difference"),
				clEnumValN(AbsVolume,"abs-volume", "Absolute volume difference"),
				clEnumValN(Hausdorff,"hausdorff", "Hausdorff distance"),
				clEnumValN(RelativeTriple, "rel-triple", "Triple including selectivity"),
				clEnumValN(AbsoluteTriple, "abs-triple", "Triple including selectivity (absolute)"),
				clEnumValEnd), cl::init(RelativeVolume));

cl::opt<unsigned> SuperNodeDepth("superdepth",
		cl::desc("Depth to descend to create aggregrated super root"),
		cl::init(0));

static void spherizeMol(OBMol& mol, vector<MolSphere>& spheres)
{
	spheres.clear();
	spheres.reserve(mol.NumAtoms());
	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
	{
		OBAtom* atom = *aitr;
		spheres.push_back(
				MolSphere(atom->x(), atom->y(), atom->z(),
						etab.GetVdwRad(atom->GetAtomicNum())));
	}

}

//do search between include and exclude
static void do_dcsearch(GSSTreeSearcher& gss, const string& includeMol,
		const string& excludeMol, const string& output, double less, double more)
{
	double dimension = gss.getDimension();
	double resolution = gss.getResolution();
	//read query molecule(s)

	//use explicit volumes
	OBConversion inmol;
	OBConversion exmol;
	inmol.SetInFormat(inmol.FormatFromExt(includeMol));
	exmol.SetInFormat(exmol.FormatFromExt(excludeMol));

	OBMol imol;
	Molecule inMol;
	if (includeMol.size() > 0 && inmol.ReadFile(&imol, includeMol))
	{
		inMol.set(imol, resolution, dimension, ProbeRadius, less);
	}
	OBMol emol;
	Molecule exMol;
	if (excludeMol.size() > 0 && exmol.ReadFile(&emol, excludeMol))
	{
		exMol.set(emol, resolution, dimension, ProbeRadius, more);
	}

	vector<Molecule> res;

	//search
	if (!ScanOnly)
		gss.dc_search(inMol, exMol, true, output.size() > 0, res);

	if (ScanCheck || ScanOnly)
	{
		vector<Molecule> res2;
		gss.dc_scan_search(inMol, exMol, true, output.size() > 0, res2);
		if (res2.size() != res.size())
		{
			cerr << "Scanning found different number: " << res2.size() << "\n";
		}
	}

	OBConversion outconv;
	outconv.SetOutFormat(outconv.FormatFromExt(output.c_str()));
	ofstream out(output.c_str());
	outconv.SetOutStream(&out);

	for (unsigned i = 0, n = res.size(); i < n; i++)
		outconv.Write(&res[i].getMol());

}

void do_nnsearch(GSSTreeSearcher& gss, const string& input, const string& output, unsigned k)
{
	//read query molecule(s)
	vector<Molecule> res;
	Molecule::iterator molitr(input, gss.getDimension(),
			gss.getResolution(), ProbeRadius);
	for (; molitr; ++molitr)
	{
		const Molecule& mol = *molitr;
		gss.nn_search(mol, k, output.size() > 0, res);
	}

	OBConversion outconv;
	outconv.SetOutFormat(outconv.FormatFromExt(output.c_str()));
	ofstream out(output.c_str());
	outconv.SetOutStream(&out);

	for (unsigned i = 0, n = res.size(); i < n; i++)
		outconv.Write(&res[i].getMol());
}

struct QInfo
{
	string str;
	CommandEnum cmd;
	Molecule in;
	Molecule ex;
	unsigned k;
	double less;
	double more;

	QInfo(): cmd(Create), k(1), less(0), more(0) {}

	QInfo(const string& s, const Molecule& i, const Molecule& e, double l, double m): str(s), cmd(DCSearch), in(i), ex(e), less(l), more(m) {}
	QInfo(const string& s, const Molecule& i, unsigned _k): str(s), cmd(NNSearch), in(i), k(_k) {}
};

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
		switch (PackerChoice)
		{
		case FullMerge:
			packer = PackerPtr(new FullMergePacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case MatchPack:
			packer = PackerPtr(
					new MatcherPacker(Pack, Knn, Sentinals, ClusterDist, QuadPack));
			break;
		case GreedyMerge:
			packer = PackerPtr(new GreedyPacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case Spectral:
			packer = PackerPtr(
					new SpectralPacker(Pack, SpectralAlg, !UseUnnorm));
			break;
		}

		setDistance(ShapeDist, MaxDimension);

		GSSLevelCreator leveler(&topdown, packer.get(), NodePack, LeafPack);

		GSSTreeCreator creator(&leveler, SuperNodeDepth);

		filesystem::path dbpath(Database.c_str());
		Molecule::iterator molitr(Input, MaxDimension, Resolution, ProbeRadius);
		if (!creator.create(dbpath, molitr, MaxDimension, Resolution))
		{
			cerr << "Error creating database\n";
			exit(1);
		}

		if (Verbose)
			creator.printStats(cout);

	}
		break;
	case NNSearch:
	{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);
		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}

		do_nnsearch(gss, Input, Output, K);
	}
		break;
	case DCSearch:
	{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}

		double dimension = gss.getDimension();
		double resolution = gss.getResolution();
		//read query molecule(s)
		if (IncludeMol.size() > 0 || ExcludeMol.size() > 0)
		{
			do_dcsearch(gss, IncludeMol, ExcludeMol, Output, LessDist, MoreDist);
		}
		else // range from single molecules
		{
			OBConversion inconv;
			inconv.SetInFormat(inconv.FormatFromExt(Input));
			ifstream in(Input.c_str());
			inconv.SetInStream(&in);

			OBMol mol;
			while (inconv.Read(&mol))
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

				vector<Molecule> res;
				//search
				if (!ScanOnly)
					gss.dc_search(
							Molecule(littlespheres, dimension, resolution),
							Molecule(bigspheres, dimension, resolution), false,
							Output.size() > 1, res);

				if (ScanCheck || ScanOnly)
				{
					vector<Molecule> res2;
					gss.dc_scan_search(
							Molecule(littlespheres, dimension, resolution),
							Molecule(bigspheres, dimension, resolution), false,
							Output.size() > 1, res2);
					if (res2.size() != res.size())
					{
						cerr << "Scanning found different number\n";
					}
				}

				OBConversion outconv;
				outconv.SetOutFormat(outconv.FormatFromExt(Output.c_str()));
				ofstream out(Output.c_str());
				outconv.SetOutStream(&out);

				for (unsigned i = 0, n = res.size(); i < n; i++)
					outconv.Write(&res[i].getMol());
			}
		}
	}
		break;
	case MolGrid:
	{
		ofstream out(Output.c_str());
		float adjust = LessDist;

		for (MolIterator mitr(Input, MaxDimension, Resolution, ProbeRadius,
				adjust); mitr; ++mitr)
		{
			MappableOctTree *tree = MappableOctTree::create(MaxDimension,
					Resolution, *mitr);
			cout << "Size of tree " << tree->bytes() << "\n";
			cout << "Volume of tree " << tree->volume() << "\n";
			cout << "Nodes of tree " << tree->nodes() << "\n";
			vector<unsigned> cnts;
			tree->countLeavesAtDepths(cnts);
			for (unsigned i = 0, n = cnts.size(); i < n; i++)
			{
				cout << i << " : " << cnts[i] << "\n";
			}
			if (out)
			{
				if (filesystem::extension(Output.c_str()) == ".raw")
					tree->dumpRawGrid(out, Resolution);
				else if (filesystem::extension(Output.c_str()) == ".mira")
					tree->dumpMiraGrid(out, Resolution);
				else
					tree->dumpGrid(out, Resolution);
			}
			free(tree);
		}
		break;
	}
	case BatchSearch:
	{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}

		//read in each line of the batch file which should be
		//cmd in_ligand in_receptor(for DC Search)
		ifstream batch(Input.c_str());
		if (!batch)
		{
			cerr << "Could not read batch file " << Input << "\n";
			exit(-1);
		}

		string line;
		while (getline(batch, line))
		{
			Timer t;
			stringstream toks(line);
			string cmd, ligand, receptor, output; //output always empty for batch
			toks >> cmd;

			if (!toks)
				break;
			if (cmd == "DCSearch")
			{
				double less = 0, more = 0;
				toks >> ligand;
				toks >> receptor;
				toks >> less;
				toks >> more;
				do_dcsearch(gss, ligand, receptor, output, less, more);
			}
			else if (cmd == "NNSearch")
			{
				unsigned k = 1;
				toks >> ligand;
				toks >> k;

				do_nnsearch(gss, ligand, output, k);
			}
			else
			{
				cerr << "Illegal command " << cmd << " in batch file.\n";
				exit(-1);
			}

			cout << "Batch " << line << " " << t.elapsed() << "\n";
		}
	}
		break;
	case BatchDB:
	{
		//setup where we run the same set of queries on a list of databases (provided in Database)
		//preproccess the queries to load in receptor/ligand structures

		//read in each line of the batch file which should be
		//cmd in_ligand in_receptor(for DC Search)
		ifstream batch(Input.c_str());
		if (!batch)
		{
			cerr << "Could not read batch file " << Input << "\n";
			exit(-1);
		}

		vector<QInfo> qinfos;
		string line;
		while (getline(batch, line))
		{
			stringstream toks(line);
			string cmd, ligand, receptor, output; //output always empty for batch
			toks >> cmd;

			if (!toks)
				break;
			if (cmd == "DCSearch")
			{
				double less = 0, more = 0;
				toks >> ligand;
				toks >> receptor;
				toks >> less;
				toks >> more;

				OBConversion inmol;
				OBConversion exmol;
				inmol.SetInFormat(inmol.FormatFromExt(ligand));
				exmol.SetInFormat(exmol.FormatFromExt(receptor));

				OBMol imol;
				Molecule inMol;
				if (ligand.size() > 0 && inmol.ReadFile(&imol, ligand))
				{
					inMol.set(imol, Resolution, MaxDimension, ProbeRadius, less);
				}

				OBMol emol;
				Molecule exMol;
				if (receptor.size() > 0 && exmol.ReadFile(&emol, receptor))
				{
					exMol.set(emol, Resolution, MaxDimension, ProbeRadius, more);
				}

				qinfos.push_back(QInfo(line, inMol, exMol, less, more));
			}
			else if (cmd == "NNSearch")
			{
				unsigned k = 1;
				toks >> ligand;
				toks >> k;

				OBConversion inmol;
				inmol.SetInFormat(inmol.FormatFromExt(ligand));

				//single mol only
				OBMol imol;
				Molecule inMol;
				if (ligand.size() > 0 && inmol.ReadFile(&imol, ligand))
				{
					inMol.set(imol, Resolution, MaxDimension, ProbeRadius, 0);
				}

				qinfos.push_back(QInfo(line, inMol, k));
			}
			else
			{
				cerr << "Illegal command " << cmd << " in batch file.\n";
				exit(-1);
			}
		}

		//read in list of databases
		ifstream dbs(Database.c_str());
		if (!dbs)
		{
			cerr << "Could not read db file list " << Database << "\n";
			exit(-1);
		}

		while (getline(dbs, line))
		{
			filesystem::path dbfile(line.c_str());
			GSSTreeSearcher gss(Verbose);

			if (!gss.load(dbfile))
			{
				cerr << "Could not read database " << line << "\n";
				exit(-1);
			}

			if(gss.getResolution() != Resolution || gss.getDimension() != MaxDimension)
			{
				cerr << "Resolution or dimension mismatch\n";
				exit(-1);
			}

			vector<Molecule> res;
			for(unsigned i = 0, n = qinfos.size(); i < n; i++)
			{
				cout << line << " " << qinfos[i].str << "\n";
				if(qinfos[i].cmd == NNSearch)
				{
					gss.nn_search(qinfos[i].in, qinfos[i].k, false, res);
				}
				else
				{
					gss.dc_search(qinfos[i].in, qinfos[i].ex, true, false, res);
				}
			}
		}

	}
		break;
	}
	return 0;
}
