/*
 * ResultMolecules.h
 *
 *  Created on: Mar 20, 2012
 *      Author: dkoes
 *
 *  Container class for reading and storing molecules from the molecular data.
 */

#ifndef RESULTMOLECULES_H_
#define RESULTMOLECULES_H_

#include <iostream>
#include <vector>
#include "PMol.h"
using namespace std;

class ResultMolecules
{
	PMolReaderMalloc reader;
	vector<PMol*> mols;

public:
	ResultMolecules()
	{

	}
	virtual ~ResultMolecules()
	{
		clear();
	}

	void clear()
	{
		for(unsigned i = 0, n = mols.size(); i < n; i++)
		{
			if(mols[i]) free(mols[i]);
		}
		mols.clear();
	}

	//add a molecule to the result set; written in pmolf format at data
	void add(const char *data)
	{
		mols.push_back((PMol*)reader.readPMol(data));
	}

	void reserve(unsigned n)
	{
		mols.reserve(n);
	}

	//number of molecules in set
	unsigned size() const
	{
		return mols.size();
	}

	//write ith mol
	void writeSDF(ostream& out, unsigned i) const
	{
		assert(i < mols.size());
		mols[i]->writeSDF(out);
	}

};

#endif /* RESULTMOLECULES_H_ */
