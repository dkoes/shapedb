/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
/*
 * MemMapped.cpp
 *
 *  Created on: Oct 19, 2011
 *      Author: dkoes
 */

#include "MemMapped.h"
#include <sys/types.h>
#include <sys/mman.h>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>

using namespace boost;

void MemMapped::clear()
{
	if(addr != NULL)
		munmap(addr, sz);
	addr = NULL;
	sz = 0;
}

bool MemMapped::map(const string& fname, bool readOnly, bool sequential, bool populate/*=false*/, bool readonce/*=false*/)
{
	if(addr != NULL)
		munmap(addr, sz);
	addr = NULL;
	unsigned flags = readOnly ? O_RDONLY : O_RDWR;
	int fd = open(fname.c_str(), flags);
	assert(fd >= 0);
#ifdef POSIX_FADV_SEQUENTIAL
	if(sequential)
		posix_fadvise(fd, 0,0,POSIX_FADV_SEQUENTIAL);
#endif
#ifdef POSIX_FADV_NOREUSE
	if(readonce)
		posix_fadvise(fd, 0, 0, POSIX_FADV_NOREUSE);
#endif
	flags = readOnly ? PROT_READ : (PROT_READ | PROT_WRITE);
	sz = filesystem::file_size(fname);
	if (sz > 0)
	{
		addr =  mmap(NULL, sz, flags, (readOnly ? MAP_PRIVATE
				: MAP_SHARED) |(populate ? MAP_POPULATE : 0) , fd, 0);
	}
	else
	{
		addr = NULL;
		sz = 0;
	}
	if ((long) addr == -1)
	{
		perror("mmap  failed: ");
		abort();
	}

	return addr != NULL;
}
