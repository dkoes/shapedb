/*
 * GSSTreeSearcher.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *      Class for loading and searching a GSSTree.
 */

#ifndef GSSTREESEARCHER_H_
#define GSSTREESEARCHER_H_

#include <boost/filesystem.hpp>
#include "GSSTreeStructures.h"
#include "Molecule.h"

using namespace boost;
using namespace std;



class GSSTreeSearcher
{
	void* objects; //memory mapped objects
	vector<void *> nodes; //memory mapped levels of tree (TODO: change to one DFO file)

public:
	GSSTreeSearcher(): objects(NULL) {}

	bool load(const filesystem::path& dbpath);
	void clear();
	~GSSTreeSearcher();
};

#endif /* GSSTREESEARCHER_H_ */
