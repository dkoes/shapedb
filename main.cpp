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
#include "molecules/Molecule.h"
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

cl::opt<CommandEnum> Command(cl::desc("Operation to perform:"), cl::Required,
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
cl::opt<PackerEnum> PackerChoice(cl::desc("Packing algorithm:"),
		cl::values(clEnumValN(FullMerge,"full-merge", "Greedy full merge."),
				clEnumValN(GreedyMerge,"greedy-merge", "Greedy iterative merge."),
				clEnumValN(MatchPack,"match-merge", "Optimal matching merging."),
				clEnumValN(Spectral, "spectral", "Spectral packing"),
				clEnumValEnd), cl::init(MatchPack));

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

cl::opt<double> LessDist("less",
		cl::desc(
				"Distance to reduce query mol by for constraint search (default 1A)."),
		cl::init(1.0));
cl::opt<double> MoreDist("more",
		cl::desc(
				"Distance to increase query mol by for constraint search (default 1A)."),
		cl::init(1.0));

cl::opt<double> ProbeRadius("probe-radius",
		cl::desc("Radius of water probe for exmol only"), cl::init(1.4));

cl::opt<double> MaxDimension("max-dim", cl::desc("Maximum dimension."),
		cl::init(64));
cl::opt<double> Resolution("resolution",
		cl::desc("Best resolution for shape database creation."), cl::init(.5));

cl::opt<string> IncludeMol("ligand",
		cl::desc("Molecule to use for minimum included volume"));
cl::opt<string> ExcludeMol("receptor",
		cl::desc("Molecule to use for excluded volume"));

cl::opt<bool> useInteractionPoints("use-interaction-points",
		cl::desc(
				"Analyze the ligand-receptor complex and generate a minimum shape centered around different interaction points."),
		cl::init(false));
cl::opt<double> interactionPointRadius("interaction-point-radius",
		cl::desc("Amount to grow interaction points into ligand shape."),
		cl::init(0));
cl::opt<double> interactionDistance("interaction-distance",
		cl::desc(
				"Distance between ligand/receptor atom centers that are considered interacting"),
		cl::init(6));
cl::opt<double> interactionMinCluster("interaction-min-cluster",
		cl::desc("Minimum size of interaction point cluster."), cl::init(3));
cl::opt<double> interactionMaxClusterDist("interaction-max-cluster-dist",
		cl::desc("Maximum span of interaction point cluster."), cl::init(4));

cl::opt<unsigned> KCenters("kcenters",
		cl::desc("number of centers for ksample-split"), cl::init(8));
cl::opt<unsigned> KSampleMult("ksamplex",
		cl::desc("multiplictive factor for ksampling"), cl::init(5));

cl::opt<unsigned> SwitchToPack("switch-to-pack",
		cl::desc("Cutoff to trigger packing in nodes and leaves"),
		cl::init(32768));

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

cl::opt<unsigned> TimeTrials("time-trials",
		cl::desc("Number of runs to get average for benchmarking"),
		cl::init(1));
cl::opt<bool> ClearCache("clear-cache",
		cl::desc("Clear file cache between each benchmarking run"),
		cl::init(false));
cl::opt<bool> ClearCacheFirst("clear-cache-first",
		cl::desc("Clear file cache before each benchmarking run"),
		cl::init(false));

cl::opt<string> SproxelColor("sproxel-color",
		cl::desc("Sproxel voxel descriptor"), cl::init("#0000ffff"));

cl::opt<bool> KeepHydrogens("h",
		cl::desc("Retain hydrogens in input molecules"), cl::init(false));

//create min and max trees from molecular data
//caller takes owner ship of tree memory
static void create_trees(GSSTreeSearcher& gss, const string& includeMol,
		const string& excludeMol,double less,
		double more, MappableOctTree*& smallTree, MappableOctTree*& bigTree)
{
	double dimension = gss.getDimension();
	double resolution = gss.getResolution();
	//read query molecule(s)

	//use explicit volumes
	Molecule::iterator imolitr(includeMol, dimension, resolution, KeepHydrogens,
			ProbeRadius);
	Molecule inMol = *imolitr;

	Molecule::iterator exmolitr(excludeMol, dimension, resolution,
			KeepHydrogens, ProbeRadius);
	Molecule exMol = *exmolitr;

	//create trees

	//resize receptor constraint
	bigTree = MappableOctTree::create(dimension, resolution,
			exMol);

	if (more > 0)
	{
		MGrid grid;
		bigTree->makeGrid(grid, resolution);
		grid.shrink(more);
		free(bigTree);
		bigTree = MappableOctTree::createFromGrid(grid);
	}
	//exmol is the receptor, so must invert
	bigTree->invert();

	//use the full shape of the ligand
	smallTree = MappableOctTree::create(dimension, resolution,
			inMol);

	//shrink if requested
	if (less > 0)
	{
		MGrid grid;
		smallTree->makeGrid(grid, resolution);
		grid.shrink(less);
		free(smallTree);
		smallTree = MappableOctTree::createFromGrid(grid);
	}

	if(useInteractionPoints)
	{
		//change small tree to just be interaction points
		MGrid igrid(MaxDimension, Resolution);
		inMol.computeInteractionGridPoints(exMol, igrid,
				interactionDistance, interactionMaxClusterDist,
				interactionMinCluster, interactionPointRadius);

		MGrid lgrid;
		smallTree->makeGrid(lgrid, Resolution);
		lgrid &= igrid;
		free(smallTree);
		smallTree = MappableOctTree::createFromGrid(lgrid);
	}

}

//do search between include and exclude
static void do_dcsearch(GSSTreeSearcher& gss, const MappableOctTree* smallTree, const MappableOctTree* bigTree, const string& output)
{
	ResultMolecules res;

	//search
	if (!ScanOnly)
		gss.dc_search(smallTree, bigTree, output.size() > 0, res);

	if (ScanCheck || ScanOnly)
	{
		ResultMolecules res2;
		gss.dc_scan_search(smallTree, bigTree, output.size() > 0,
				res2);
		if (res2.size() != res.size())
		{
			cerr << "Scanning found different number: " << res2.size() << "\n";
		}
	}

	if (output.size() > 0)
	{
		ofstream outmols(output.c_str());
		for (unsigned i = 0, n = res.size(); i < n; i++)
			res.writeSDF(outmols, i);
	}

}

void do_nnsearch(GSSTreeSearcher& gss, const string& input,
		const string& output, unsigned k)
{
	//read query molecule(s)
	ResultMolecules res;
	Molecule::iterator molitr(input, gss.getDimension(), gss.getResolution(),
			KeepHydrogens, ProbeRadius);
	for (; molitr; ++molitr)
	{
		const Molecule& mol = *molitr;
		if(k < 1) //scan
		{
			gss.nn_scan(mol, output.size() > 0, res);
		}
		else
		{
			gss.nn_search(mol, k, output.size() > 0, res);
		}
	}

	ofstream out(output.c_str());
	for (unsigned i = 0, n = res.size(); i < n; i++)
		res.writeSDF(out, i);
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

	QInfo() :
			cmd(Create), k(1), less(0), more(0)
	{
	}

	QInfo(const string& s, const Molecule& i, const Molecule& e, double l,
			double m) :
			str(s), cmd(DCSearch), in(i), ex(e), less(l), more(m)
	{
	}
	QInfo(const string& s, const Molecule& i, unsigned _k) :
			str(s), cmd(NNSearch), in(i), k(_k)
	{
	}
};

static void dumpTree(MappableOctTree *tree, ostream& out)
{

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
		else if (filesystem::extension(Output.c_str()) == ".csv")
			tree->dumpSproxelGrid(out, Resolution, SproxelColor);
		else
			tree->dumpGrid(out, Resolution);
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
		KSamplePartitioner topdown(KCenters, KSampleMult,
				KSamplePartitioner::AveCenter, SwitchToPack);

		PackerPtr packer;
		switch (PackerChoice)
		{
		case FullMerge:
			packer = PackerPtr(
					new FullMergePacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case MatchPack:
			packer = PackerPtr(
					new MatcherPacker(Pack, Knn, Sentinals, ClusterDist));
			break;
		case GreedyMerge:
			packer = PackerPtr(
					new GreedyPacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case Spectral:
			packer = PackerPtr(
					new SpectralPacker(Pack, SpectralAlg, !UseUnnorm));
			break;
		}

		setDistance(ShapeDist, MaxDimension);

		GSSLevelCreator leveler(&topdown, packer.get(), SwitchToPack,
				SwitchToPack);

		GSSTreeCreator creator(&leveler, SuperNodeDepth);

		filesystem::path dbpath(Database.c_str());
		Molecule::iterator molitr(Input, MaxDimension, Resolution,
				KeepHydrogens, ProbeRadius);
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

		if(ExcludeMol.size() == 0)
		{
			do_nnsearch(gss, IncludeMol, Output, K);
		}
		else
		{
			ResultMolecules res;
			MappableOctTree *smallTree = NULL, *bigTree = NULL;
			create_trees(gss, IncludeMol, ExcludeMol, LessDist, MoreDist, smallTree, bigTree);
			if(K < 1) //scan
			{
				gss.nn_scan(smallTree, bigTree, Output.size() > 0, res);
			}
			else
			{
				gss.nn_search(smallTree, bigTree, K, Output.size() > 0, res);
			}
			if(Output.size() > 0)
			{
				ofstream out(Output.c_str());
				for (unsigned i = 0, n = res.size(); i < n; i++)
					res.writeSDF(out, i);
			}
		}
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
			MappableOctTree *smallTree = NULL, *bigTree = NULL;
			create_trees(gss, IncludeMol, ExcludeMol, LessDist, MoreDist, smallTree, bigTree);
			do_dcsearch(gss, smallTree, bigTree, Output);
			free(smallTree);
			free(bigTree);
		}
		else // range from single molecules
		{
			for (Molecule::iterator inmols(Input, dimension, resolution,
					KeepHydrogens, ProbeRadius); inmols; ++inmols)
			{
				ResultMolecules res;
				//search
				MappableOctTree* tree = MappableOctTree::create(dimension, resolution, *inmols);
				MGrid lgrid(dimension,resolution);
				tree->makeGrid(lgrid, resolution);
				free(tree);

				MGrid lgrid2(lgrid);

				//resize
				lgrid.shrink(LessDist);
				lgrid2.grow(MoreDist);

				//create bounding trees
				MappableOctTree *smallTree = MappableOctTree::createFromGrid(lgrid);
				MappableOctTree *bigTree = MappableOctTree::createFromGrid(lgrid2);

				if (!ScanOnly)
					gss.dc_search(smallTree, bigTree, Output.size() > 1, res);

				if (ScanCheck || ScanOnly)
				{
					ResultMolecules res2;
					gss.dc_scan_search(smallTree,bigTree, Output.size() > 1, res2);
					if (res2.size() != res.size())
					{
						cerr << "Scanning found different number\n";
					}
				}

				free(smallTree);
				free(bigTree);

				ofstream outmols(Output.c_str());
				for (unsigned i = 0, n = res.size(); i < n; i++)
					res.writeSDF(outmols, i);
			}
		}
	}
		break;
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
			stringstream toks(line);
			string cmd, ligand, receptor, output; //output always empty for batch
			toks >> cmd;

			if (!toks)
				break;

			double less = 0, more = 0;
			unsigned k = 1;

			if (cmd == "DCSearch")
			{
				toks >> ligand;
				toks >> receptor;
				toks >> less;
				toks >> more;
			}
			else if (cmd == "NNSearch")
			{
				toks >> ligand;
				toks >> k;
			}

			vector<double> times;

			if (ClearCacheFirst)
				std::system("clearfilecache");

			for (unsigned i = 0; i < TimeTrials; i++)
			{
				if (ClearCache)
				{
					std::system("clearfilecache");
				}

				Timer t;

				if (cmd == "DCSearch")
				{
					MappableOctTree *smallTree = NULL, *bigTree = NULL;
					create_trees(gss, ligand, receptor, less, more, smallTree, bigTree);
					do_dcsearch(gss, smallTree, bigTree, output);
					free(smallTree); free(bigTree);
				}
				else if (cmd == "NNSearch")
				{
					do_nnsearch(gss, ligand, output, k);
				}
				else
				{
					cerr << "Illegal command " << cmd << " in batch file.\n";
					exit(-1);
				}

				times.push_back(t.elapsed());
			}

			cout << "Batch " << line;
			for (unsigned i = 0, n = times.size(); i < n; i++)
			{
				cout << " " << times[i];
			}
			cout << "\n";
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

				Molecule::iterator inmol(ligand, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule inMol = *inmol;

				Molecule::iterator exmol(receptor, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule exMol = *exmol;

				qinfos.push_back(QInfo(line, inMol, exMol, less, more));
			}
			else if (cmd == "NNSearch")
			{
				unsigned k = 1;
				toks >> ligand;
				toks >> k;

				Molecule::iterator inmol(ligand, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule inMol = *inmol;

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

			if (gss.getResolution() != Resolution
					|| gss.getDimension() != MaxDimension)
			{
				cerr << "Resolution or dimension mismatch\n";
				exit(-1);
			}

			ResultMolecules res;
			for (unsigned i = 0, n = qinfos.size(); i < n; i++)
			{
				cout << line << " " << qinfos[i].str << "\n";
				if (qinfos[i].cmd == NNSearch)
				{
					gss.nn_search(qinfos[i].in, qinfos[i].k, false, res);
				}
				else
				{
					gss.dc_search(qinfos[i].in, qinfos[i].ex, qinfos[i].less,
							qinfos[i].more, true, false, res);
				}
			}
		}

	}
		break;
	case MolGrid:
	{
		ofstream out(Output.c_str());

		if (useInteractionPoints)
		{
			if (IncludeMol.size() == 0 || ExcludeMol.size() == 0)
			{
				cerr
						<< "Need both ligand and receptor for interaction points.\n";
				exit(-1);
			}
			Molecule::iterator imolitr(IncludeMol, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius);
			Molecule ligand = *imolitr;

			Molecule::iterator exmolitr(ExcludeMol, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius);
			Molecule receptor = *exmolitr;

			MGrid igrid(MaxDimension, Resolution);
			ligand.computeInteractionGridPoints(receptor, igrid,
					interactionDistance, interactionMaxClusterDist,
					interactionMinCluster, interactionPointRadius);

			MappableOctTree *tree = MappableOctTree::create(MaxDimension,
					Resolution, ligand);
			MGrid lgrid;
			tree->makeGrid(lgrid, Resolution);
			lgrid.shrink(LessDist);

			lgrid &= igrid;

			free(tree);
			tree = MappableOctTree::createFromGrid(lgrid);

			dumpTree(tree, out);
			free(tree);
		}
		else
		{
			for (Molecule::iterator mitr(Input, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius); mitr; ++mitr)
			{
				MappableOctTree *tree = MappableOctTree::create(MaxDimension,
						Resolution, *mitr);

				if (LessDist >= 0)
				{
					MGrid grid;
					tree->makeGrid(grid, Resolution);

					grid.shrink(LessDist);
					free(tree);
					tree = MappableOctTree::createFromGrid(grid);
				}

				dumpTree(tree, out);
				free(tree);
			}
		}
		break;
	}
	}
	return 0;
}
