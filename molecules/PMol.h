/*
 * PMol.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 *
 *  This is specialized class for storing small molecule information.
 *  It has a fairly compressed binary output.  There is a single conformer
 *  per each molecule since this makes lookup faster.
 */

#ifndef PMOL_H_
#define PMOL_H_


#include <boost/unordered_map.hpp>
#include <openbabel/mol.h>
#include <vector>
#include <rdkit/GraphMol/ROMol.h>
#include "../Orienter.h"

using namespace std;

#define MAX_BONDS 3


struct FloatCoord
{
	float x;
	float y;
	float z;

	FloatCoord(float a, float b, float c) :
		x(a), y(b), z(c)
	{
	}
	FloatCoord() :
		x(0), y(0), z(0)
	{
	}

};
struct Property
{
	unsigned char atom;
	char value;

};

/* Class for creating an PMol (from a more expressive mol) and writing out */
class PMolCreator
{
	//atoms are stored grouped by typ
	struct AtomGroup
	{
		unsigned short startIndex; //for convenience, not output to file
		unsigned char atomic_number;
		vector<FloatCoord> coords;

		AtomGroup(unsigned anum) :
			startIndex(0), atomic_number(anum)
		{
		}
	};

	//Bounds are grouped by type (single, double, triple)
	//and stored as adjacency lists

	//adjacencly lists, segregated by bond type
	//directly indexed by atom index; one directional
	vector<vector<unsigned char> > bonds[MAX_BONDS];
	vector<AtomGroup> atoms;

	struct Property
	{
		unsigned char atom;
		char value;

		Property(unsigned ai, char v) :
			atom(ai), value(v)
		{
		}
	};

	//additional properties
	vector<Property> iso;
	vector<Property> chg;

	string name;
	unsigned numAtoms;
	unsigned nSrcs;
	unsigned bndSize[MAX_BONDS];
	unsigned nDsts;
public:
	PMolCreator()
	{
	}
	PMolCreator(OpenBabel::OBMol& mol, bool deleteH=false) :
		numAtoms(0), nSrcs(0), nDsts(0)
	{
		memset(bndSize, 0, sizeof(bndSize));
		copyFrom(mol, deleteH);
	}
	PMolCreator(RDKit::ROMol& mol, bool deleteH=false) :
		numAtoms(0), nSrcs(0), nDsts(0)
	{
		memset(bndSize, 0, sizeof(bndSize));
		copyFrom(mol, deleteH);
	}
	~PMolCreator()
	{
	}

	void copyFrom(OpenBabel::OBMol& mol, bool deleteH=false);
	void copyFrom(RDKit::ROMol& mol, bool deleteH=false);

	//write custom binary data
	bool writeBinary(ostream& out) const;

};

struct AdjList
{
	unsigned char src;
	unsigned char nDsts;
	unsigned char dsts[];
};

#define PMOLHEADER_MAX (256)
struct PMolHeader
{
	unsigned char nAtoms;
	unsigned char nAtomTypes;
	unsigned char nISO;
	unsigned char nCHG;
	unsigned char nBnds; //this could actually be computed, so is available space
	unsigned char adjListSize[MAX_BONDS];
	FloatCoord coords[];
};

struct AtomTypeCnts
{
	unsigned char atomic_number;
	unsigned char cnt;
};

struct ASDDataItem
{
	string tag;
	string value;

	ASDDataItem(const string& t, const string& v) :
		tag(t), value(v)
	{
	}
};

//Created by PMolReader, can output to regular mol formats, as well as be reoriented
//Data is stored of the end of the structure and we store pointers into this data
class PMol
{
	AtomTypeCnts *atomtypes; //nAtomTypes
	Property *iso; //nISO many
	Property *chg; //nCHG many
	char *adjlists[MAX_BONDS]; //cast to adjlist and iterate carefully
	char *name;
	PMolHeader header;
	char buffer[];

	PMol() :
		atomtypes(NULL), iso(NULL), chg(NULL)
	{
		memset(adjlists, 0, sizeof(adjlists));
	}
	friend class PMolReader;

	//assume data is already copied in, setup pointers
	//return final offset
	unsigned setup();
public:

	const char *getTitle() const
	{
		return name;
	}
	double getMolWeight() const;

	//write sdf with associated meta data
	//rotate/translate points if requested
	void writeSDF(ostream& out, const vector<ASDDataItem>& sddata,
			const Orienter& rms);

	void writeSDF(ostream& out, const Orienter& rms)
	{
    vector<ASDDataItem> dummy;
		writeSDF(out, dummy, rms);
	}

  void writeSDF(ostream& out, const vector<ASDDataItem>& sddata)
  {
    Orienter rdummy;
    writeSDF(out, sddata, rdummy);
  }

	void writeSDF(ostream& out)
	{
		vector<ASDDataItem> dummy;
		Orienter rdummy;
		writeSDF(out, dummy, rdummy);
	}

	//mutating rotation/translation
	void rotate(const double *rotation);
	void translate(const double *translation);
};

//generated PMol's
class PMolReader
{
protected:
	//read into an already allocated buffer
	virtual void* allocate(unsigned size) = 0;
public:
	virtual ~PMolReader() {}
	virtual PMol* readPMol(const char *data); //return pmol at data
	virtual PMol* readPMol(FILE *f); //return amol at current pos

};


//malloc's the required memory, must be freed by caller
class PMolReaderMalloc: public PMolReader
{

protected:
	virtual void* allocate(unsigned size);

};

//owns and maintains a buffer suitable for a single mol
//each read will overwrite previous - do NOT share between threads
class PMolReaderSingleAlloc: public PMolReader
{
	void *buffer;
	unsigned bsize;
protected:
	virtual void* allocate(unsigned size);
public:
	PMolReaderSingleAlloc():bsize(2048)
	{
		buffer = malloc(bsize);
	}

	~PMolReaderSingleAlloc()
	{
		free(buffer);
	}

};


//uses a fixed size buffer (so can be stack allocated)
//dies if buffer isn't big enough
class PMolReaderSimple: public PMolReader
{
	static const int buffsize = 8192;
	char buffer[buffsize];
protected:
	virtual void* allocate(unsigned size) {
		assert(size < buffsize);
		return (void*) buffer;
	}
public:
	PMolReaderSimple()
	{
	}

	~PMolReaderSimple()
	{
	}

};

#endif /* PMOL_H_ */
