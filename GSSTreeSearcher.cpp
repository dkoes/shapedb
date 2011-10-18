/*
 * GSSTreeSearcher.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#include "GSSTreeSearcher.h"

//load a gsstree database by mmapping files, return true if successfull
bool GSSTreeSearcher::load(const filesystem::path& dbpath)
{
	clear();
	//memory map files form db directory


	return true;
}

void GSSTreeSearcher::clear()
{
	//unmap if necessary
	abort();
}

GSSTreeSearcher::~GSSTreeSearcher()
{
	clear();
}

