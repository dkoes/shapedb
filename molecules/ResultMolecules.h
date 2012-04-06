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
#include <boost/lexical_cast.hpp>
#include "PMol.h"
using namespace std;

class ResultMolecules
{
	PMolReaderMalloc reader;
	vector<PMol*> mols;
	vector<double> scores;
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
	//support a single sddata entry for the goodness of the result
	void add(const char *data, double score)
	{
		mols.push_back((PMol*)reader.readPMol(data));
		scores.push_back(score);
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
		vector<ASDDataItem> data;
		data.push_back(ASDDataItem("score",boost::lexical_cast<string>(scores[i])));
		mols[i]->writeSDF(out, data);
	}

};

#endif /* RESULTMOLECULES_H_ */
