/**
 * Utility program to read data output
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "util.hpp"
#include "ser_io.h"

using namespace std;

int main( int argc, char** argv )
{

	int nBlock = 0;
	unsigned int* blockPointerList = NULL;
	char *fileName = argv[1];
	 
	 // Open file
	FILE* dataFile = fopen( fileName, "rb" );
	assert( dataFile != NULL );

	// Initialize IO class
	SER_IO<float> *ser_io = new SER_IO<float>( 3, false );

	// Read footer
	ser_io->ReadFooter( dataFile, blockPointerList, nBlock );
	printf( "The file has %d blocks\n", nBlock );

	// Allocate memory for entropy list
	float globalEntropyList[nBlock];

	//Read blocks
	ser_io->ReadBlocks( dataFile, blockPointerList, &globalEntropyList[0], nBlock );
	for( int iB=0; iB<nBlock; iB++ )
		printf( "%d: Block offset: %d Entopy: %f\n", iB, blockPointerList[iB], globalEntropyList[iB] );


	// Close file
	fclose( dataFile );

	return 0;

	
 }
