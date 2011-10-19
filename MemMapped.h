/*
 * MemMapped.h
 *
 *  Created on: Oct 19, 2011
 *      Author: dkoes
 *
 *  A wrapper for a memory mapped file.  I use this instead of boost
 *  since this way I have more control over the posix settings and
 *  a mroe convenient interface
 */

#ifndef MEMMAPPED_H_
#define MEMMAPPED_H_

#include "GSSTypes.h"
#include <string>

class MemMapped
{
	void *addr;
	unsigned long size;
public:
	MemMapped(): addr(NULL), size(0) {}
	~MemMapped() {} //does not unmap, must explicitly clear

	//map a file into memory
	bool map(const string& fname, bool readOnly, bool sequential, bool populate=false, bool readonce=false);

	const unsigned char& operator[](unsigned i) const { return ((unsigned char*)addr)[i]; }
	operator void*() const { return addr; }

	const char * begin() const { return (const char*)addr; }
	void clear();
};

#endif /* MEMMAPPED_H_ */
