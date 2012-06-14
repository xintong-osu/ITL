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

	int nBlock = 1024;
	unsigned int* blockPointerList = NULL;
	char *fileName = argv[1];
	char *outFileName = argv[2];
	//int nX = 64;
	//int nY = 64;
	//int nZ = 64;
	int nX = 126;
	int nY = 126;
	int nZ = 512;

	// Allocate memory for full entropy field
	float *entropyField = new float[nX*nY*nZ];
	 
	 // Open file
	FILE* dataFile = fopen( fileName, "rb" );
	assert( dataFile != NULL );

	// Initialize IO class
	SER_IO<float> *ser_io = new SER_IO<float>( 3, false );

	// Read footer
	ser_io->ReadFooter( dataFile, blockPointerList, nBlock );
	printf( "The file has %d blocks\n", nBlock );

	// Allocate memory for block Limits
	float **blockLimitList = new float*[nBlock];
	for( int iB=0; iB<nBlock; iB++ )
		blockLimitList[iB] = new float[6];

	//Read block limits
	ser_io->ReadBlockLimits( dataFile, blockPointerList, blockLimitList, nBlock );
	for( int iB=0; iB<nBlock; iB++ )
		printf( "%d: Block offset: %d Limit: %f %f %f %f %f %f\n", iB, blockPointerList[iB],
				blockLimitList[iB][0], blockLimitList[iB][1], blockLimitList[iB][2],
				blockLimitList[iB][3], blockLimitList[iB][4], blockLimitList[iB][5] );

	// Allocate memory for entropy fields for each block
	float *entropyFieldList[nBlock];
	int blockSizeList[nBlock];
	for( int iB=0; iB<nBlock; iB++ )
	{
		blockSizeList[iB] = (int)( 	(blockLimitList[iB][3] - blockLimitList[iB][0] + 1) *
									(blockLimitList[iB][4] - blockLimitList[iB][1] + 1) *
									(blockLimitList[iB][5] - blockLimitList[iB][2] + 1) );
		entropyFieldList[iB] = new float[blockSizeList[iB]];
	}

	//Read blocks
	ser_io->ReadBlocks( dataFile, blockPointerList, &entropyFieldList[0], blockSizeList, nBlock );
	printf( "All blocks read\n" );

	// Close file
	fclose( dataFile );

	// Combine all blocks to a global field
	int globalIndex = 0;
	float min = 100000;
	float max = -10000000;
	for( int iB=0; iB<nBlock; iB++ )
	{
		int localIndex = 0;
		for( int z=blockLimitList[iB][2]; z<=blockLimitList[iB][5]; z++ )
		{
			for( int y=blockLimitList[iB][1]; y<=blockLimitList[iB][4]; y++ )
			{
				for( int x=blockLimitList[iB][0]; x<=blockLimitList[iB][3]; x++ )
				{
					//globalIndex = z*nX*nY + x*nY + y;
					globalIndex = z*nX*nY + y*nX + x;
					entropyField[globalIndex] = entropyFieldList[iB][localIndex];

					if( entropyField[globalIndex] < min ) min = entropyField[globalIndex];
					if( entropyField[globalIndex] > max ) max = entropyField[globalIndex];

					localIndex ++;
				}
			}
		}
	}

	//printf( "%f %f\n", min, max );

	// Open file for writing
	FILE* outFile = fopen( outFileName, "wb" );
	assert( outFile != NULL );

	fwrite( &nX, sizeof(int), 1, outFile );
	fwrite( &nY, sizeof(int), 1, outFile );
	fwrite( &nZ, sizeof(int), 1, outFile );
	fwrite( entropyField, sizeof(float), nX*nY*nZ, outFile );

	fclose( outFile );



	return 0;

	
 }
