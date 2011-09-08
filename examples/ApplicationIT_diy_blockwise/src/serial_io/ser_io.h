//------------------------------------------------------------------------------
//
// Serial io class for reading blockwise global entropy values
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
// 
// All rights reserved. May not be used, modified, or copied
// without permission
//
//--------------------------------------------------------------------------

#ifndef _SER_IO
#define _SER_IO

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include "util.hpp"
#include "zlib.h"

# define HDR_ELEMENTS 0
# define HDR_SIZE 0

//------------------------------------------------------------------------------
//
// definition
//
//----------------------------------------------------------------------------
template<typename T> 
class SER_IO
{

public:

	SER_IO( int dim, bool swap_bytes = false )
	{ 
		this->dim = dim;
		this->swap_bytes = swap_bytes;
  	}

	~SER_IO(){};

	void ReadFooter( FILE*& fd, unsigned int*& ftr, int& tb );
	void ReadHeader( FILE *fd, unsigned int *hdr, size_t ofst );
	void ReadBlocks( FILE *fd, unsigned int *blockPtrList, T *globalEntropyList, int nBlock );
	void ReadBlockData( FILE *fd, unsigned int blockPtr, T* blockEntropy );

	int dim; 		// number of dimensions in the dataset
	bool swap_bytes; 	// whether to swap bytes for endian conversion

};

//------------------------------------------------------------------------------
//
// implementation
//
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//
// reads the file footer
// footer in file is always ordered by global block id
// output footer is in the same order
//
// fd: open file
// ftr: footer data (output)
// tb: total number of blocks in the file (output)
//
// side effects: allocates ftr
//
template<typename T> 
void SER_IO<T>::ReadFooter( FILE*& fd, unsigned int*& ftr, int& tb )
{
	int i;

	// Move to the end of the file to read the number of blocks
	fseek( fd, -sizeof(uint32_t), SEEK_END );

	// Read total number of blocks	
	assert( fread(&tb, sizeof(uint32_t), 1, fd) == 1 );
	if ( swap_bytes )
		swap( (char *)&tb, 1, sizeof(uint32_t) );

	// Read pointers (local offsets) to each block
	if (tb > 0)
 	{
		ftr = new uint32_t[tb];
		fseek(fd, -(tb + 1) * sizeof(uint32_t), SEEK_END);
		assert(fread(ftr, sizeof(uint32_t), tb, fd) == tb);
		if ( swap_bytes )
			swap( (char *)ftr, tb, sizeof(uint32_t) );
	}

}
//----------------------------------------------------------------------------
//
// reads the header for one block from a file
//
// fd: open file
// hdr: allocated header data
// ofst: location in file of the header (bytes)
//
template<typename T> 
void SER_IO<T>::ReadHeader( FILE *fd, unsigned int *hdr, size_t ofst )
{
	// Move to the start of the block
	fseek( fd, ofst, SEEK_SET );

	// Read block header. if any
	assert(fread( hdr, sizeof(unsigned int), HDR_ELEMENTS, fd) == HDR_ELEMENTS );
	if (swap_bytes)
		swap((char *)hdr, HDR_ELEMENTS, sizeof(unsigned int));

}
//----------------------------------------------------------------------------
//
// Reads data iteratively from the blocks from a file
//
// fd: open file
// hdr: allocated header data
// ofst: location in file of the header (bytes)
//
template<typename T>
void SER_IO<T>::ReadBlocks( FILE *fd, unsigned int *blockPtrList, T *globalEntropyList, int nBlock )
{
	for( int iB = 0; iB < nBlock; iB++ )
	{
		ReadBlockData( fd, blockPtrList[iB], globalEntropyList+iB );
	}// end for
}
//----------------------------------------------------------------------------
//
// reads data for one block from a file
//
// fd: open file
// blockPtr: Pointer to the block in file
//
template<typename T>
void SER_IO<T>::ReadBlockData( FILE *fd, unsigned int blockPtr, T *blockEntropy )
{
	// Move to the start of the block and skip the header, if any
	fseek( fd, (long int)(blockPtr + HDR_SIZE), SEEK_SET );

	// Read block data
	assert( fread( blockEntropy, sizeof(T), 1, fd ) == 1 );
	if (swap_bytes)
		swap( (char *)blockEntropy, 1, sizeof(T) );
}

#endif
