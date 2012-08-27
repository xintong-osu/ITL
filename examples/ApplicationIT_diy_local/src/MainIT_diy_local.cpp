/**
 * @file MainIT_diy_local.cpp
 * Application program for local entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon 
 */

//#define BYTE_SWAP

#include <mpi.h>

#include "diy.h"
#include "util.hpp"
#include "assignment.hpp"
#include "blocking.hpp"
#include "io.hpp"
#include "merge.hpp"

#include "bil.h"

#include "ITL_header.h"
#include "ITL_base.h"
#include "ITL_util.h"
#include "ITL_ioutil.h"
#include "ITL_vectormatrix.h"
#include "ITL_histogramconstants.h"
#include "ITL_histogram.h"
#include "ITL_histogrammapper.h"
#include "ITL_field_regular.h"
#include "ITL_globalentropy.h"
#include "ITL_localentropy.h"

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

char* inputFieldFile = NULL; 	
char* outFile = NULL;
char* patchFile = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR>* histMapper_scalar = NULL;
ITL_histogrammapper<VECTOR3>* histMapper_vector = NULL;

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_field_regular<int> *binField = NULL;

ITL_globalentropy<SCALAR> *globalEntropyComputer_scalar = NULL;
ITL_localentropy<SCALAR> *localEntropyComputer_scalar = NULL;
ITL_localentropy<VECTOR3> *localEntropyComputer_vector = NULL;

int dataSize[3];
int blockSize[3];
int tot_blocks = 512;

double execTime[3];
clock_t starttime, endtime;

//
// user-defined callback function for creating an MPI datatype for the
//   received item
//
// item: pointer to the item
// char * is used as a generic pointers to bytes, not necessarily to strings
// abs_addr: whether offsets in the MPI datatype are absolutely addressed via
//  MPI_BOTTOM, or relatively addressed via the start of the item
//  true: uses absolute addressing, false: uses relative addressing
//
// side effects: creates & commits the MPI datatype
//
// returns: pointer to the datatype
//
void*
CreateWriteType( void *item, int lid, DIY_Datatype *dtype )
{
	int min[3], size[3]; // block extents
	DIY_Block_starts_sizes( lid, min, size );
	int block_size = size[0] * size[1] * size[2];
	fprintf( stderr, "Block id: %d and size: %d\n", lid, block_size );
	DIY_Create_vector_datatype( block_size + 6, 1, DIY_FLOAT, dtype );
	return item;
}

/**
 * Main function.
 * Program starts from here.
 */
int main( int argc, char** argv )
{
	int numProcs;
	int rank;

	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	int nDim;
	int nBin;
	
	int given[3] = {0, 0, 0};

	int fieldType = -1;
	int method = 0;
	int verboseMode = 1;

	float lowF[3];
	float highF[3];
	int lowPad[3] = {0, 0, 0};
	int highPad[3] = {0, 0, 0};
	int paddedLow[3] = {0, 0, 0};
	int paddedHigh[3] = {0, 0, 0};
	int neighborhoodSizeArray[3] = {0, 0, 0};
	int neighborhoodSize = 0;
	float histogramLowEnd = 0;
	float histogramHighEnd = 0;

	int nblocks;						// My local number of blocks

	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	inputFieldFile = ITL_util<float>::getArgWithName( "inputField", &argNames, &argValues );
	outFile =  ITL_util<float>::getArgWithName( "outFile", &argNames, &argValues );
	patchFile =  ITL_util<float>::getArgWithName( "patchFile", &argNames, &argValues );
	fieldType = atoi( ITL_util<float>::getArgWithName( "fieldType", &argNames, &argValues ) );
	nDim = atoi( ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	nBin = atoi( ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
	dataSize[0] = atoi( ITL_util<float>::getArgWithName( "nX", &argNames, &argValues ) );
	dataSize[1] = atoi( ITL_util<float>::getArgWithName( "nY", &argNames, &argValues ) );
	dataSize[2] = atoi( ITL_util<float>::getArgWithName( "nZ", &argNames, &argValues ) );
	tot_blocks = atoi( ITL_util<float>::getArgWithName( "nBlock", &argNames, &argValues ) );
	neighborhoodSize = atoi( ITL_util<float>::getArgWithName( "neighborhoodSize", &argNames, &argValues ) );
	histogramLowEnd = atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramHighEnd = atof( ITL_util<float>::getArgWithName( "histHigh", &argNames, &argValues ) );
	verboseMode = atoi( ITL_util<float>::getArgWithName( "verbose", &argNames, &argValues ) );

	// Initialize ITL
	ITL_base::ITL_init();

	// Initialize histogram
	histogram = new ITL_histogram( patchFile, nBin );

	// Initialize data to histogram converter
	if( fieldType == 0 )
		histMapper_scalar = new ITL_histogrammapper<SCALAR>( histogram );
	else if( fieldType == 1 )
		histMapper_vector = new ITL_histogrammapper<VECTOR3>( histogram );

	// Allocate memory for pointers that will hold block data
	MPI_Datatype complex;
	if( fieldType == 1 )
	{
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  		MPI_Type_commit( &complex );
	}

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Decompose domain (first without ghost layers)
	DIY_Decompose( 0, 0, 0, given );

	// Serially visit the blocks (?)
	int* diy_min_noghost = new int[3*nblocks];
	int* diy_max_noghost = new int[3*nblocks];
	int* diy_size_noghost = new int[3*nblocks];

	starttime = ITL_util<float>::startTimer();
	for (int i = 0; i < nblocks; i++)
	{
		DIY_Block_starts_sizes(i, &diy_min_noghost[3*i], &diy_size_noghost[3*i] );

		// print the block bounds
		for (int j = 0; j < 3; j++)
			diy_max_noghost[3*i+j] = diy_min_noghost[3*i+j] + diy_size_noghost[3*i+j] - 1;

		if( verboseMode == 1 )
			printf("process rank = %d "
				"block local id = %d "
				"min = [%d %d %d] "
				"max = [%d %d %d] "
				"size = [%d %d %d]\n",
				rank, i,
				diy_min_noghost[3*i], diy_min_noghost[3*i+1], diy_min_noghost[3*i+2],
				diy_max_noghost[3*i], diy_max_noghost[3*i+1], diy_max_noghost[3*i+2],
				diy_size_noghost[3*i], diy_size_noghost[3*i+1], diy_size_noghost[3*i+2] );
	}

	// Decompose domain (with ghost layers) 
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	DIY_Decompose( 0, neighborhoodSize, 0, given );

	// Allocate memory for pointers that will hold block data
	float* data[nblocks];	
	VECTOR3* vectordata[nblocks];		
	if( fieldType == 0 ) memset( data, 0, sizeof(float*) * nblocks );
	if( fieldType == 1 ) memset( vectordata, 0, sizeof(VECTOR3*) * nblocks );

	// Serially visit the blocks (?)
	int* diy_min = new int[3*nblocks];
	int* diy_max = new int[3*nblocks];
	int* diy_size = new int[3*nblocks];

	starttime = ITL_util<float>::startTimer();	
	for (int i = 0; i < nblocks; i++)
	{ 
		DIY_Block_starts_sizes(i, &diy_min[3*i], &diy_size[3*i] );
	
		// post a read for the block
		if( fieldType == 0 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, DIY_FLOAT, (void**)&(data[i]));
		if( fieldType == 1 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, complex, (void**)&(vectordata[i]));

		// print the block bounds
		for (int j = 0; j < 3; j++)
			diy_max[3*i+j] = diy_min[3*i+j] + diy_size[3*i+j] - 1;

		if( verboseMode == 1 )
			printf("process rank = %d "
				"block local id = %d " 
				"min = [%d %d %d] "
				"max = [%d %d %d] "
				"size = [%d %d %d]\n", 
				rank, i, 
				diy_min[3*i], diy_min[3*i+1], diy_min[3*i+2], 
				diy_max[3*i], diy_max[3*i+1], diy_max[3*i+2], 
				diy_size[3*i], diy_size[3*i+1], diy_size[3*i+2] );
	}

	// Read actual data (everyone synchronizes after reading data)
	DIY_Read_data_all();
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	#ifdef BYTE_SWAP
	int nBlockElem = 0;
	// For byte swapping, also need to do this for ghost layers
	for( int k=0; k<nblocks; k++ )
	{
		if( fieldType == 0 )
		{
			nBlockElem = diy_size[3*k]*diy_size[3*k+1]*diy_size[3*k+2];
			swap((char *)&data[k][0], nBlockElem, sizeof( SCALAR ));

			SCALAR m = ITL_util<SCALAR>::Min( data[k], nBlockElem );
			SCALAR M = ITL_util<SCALAR>::Max( data[k], nBlockElem );
			cout << nBlockElem << " " << m << " " << M << endl;
		}
		else
		{
			nBlockElem = diy_size[3*k]*diy_size[3*k+1]*diy_size[3*k+2];
			swap((char *)&vectordata[k][0], nBlockElem*3, sizeof( SCALAR ));
		}
	}
	#endif

	// Allocate memory for storing local entropy fields for each block 		
	SCALAR *localEntropyList[nblocks];

	// Scan through the blocks, compute local entropy fields 
	starttime = ITL_util<float>::startTimer();	
	int nBlockElem = 0;
  	int *freqList[nblocks];
	for( int k=0; k<nblocks; k++ )
	{

		//lowF[0] = 0; 	highF[0] = diy_max[3*k] - diy_min[3*k];
		//lowF[1] = 0; 	highF[1] = diy_max[3*k+1] - diy_min[3*k+1];
		//lowF[2] = 0; 	highF[2] = diy_max[3*k+2] - diy_min[3*k+2];
		lowF[0] = diy_min_noghost[3*k]; 		highF[0] = diy_max_noghost[3*k];
		lowF[1] = diy_min_noghost[3*k+1]; 		highF[1] = diy_max_noghost[3*k+1];
		lowF[2] = diy_min_noghost[3*k+2]; 		highF[2] = diy_max_noghost[3*k+2];
		//dataSize[0] = (int)( highF[0] - lowF[0] + 1 );
		//dataSize[1] = (int)( highF[1] - lowF[1] + 1 );
		//dataSize[2] = (int)( highF[2] - lowF[2] + 1 );

		fprintf( stderr, "Block size (no ghost): %d %d %d\n", dataSize[0], dataSize[1], dataSize[2] );
		fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																highF[0], highF[1], highF[2] );

		if( fieldType == 0 )
		{
			// Initialize ITL scalar field with block data
			//cout << "0" << endl;
			for( int i=0; i<nDim; i++ )
			{
				neighborhoodSizeArray[i] = neighborhoodSize;
				paddedLow[i] = ITL_util<int>::clamp( (int)lowF[i] - neighborhoodSizeArray[i], 0, dataSize[i]-1 );
				paddedHigh[i] = ITL_util<int>::clamp( (int)highF[i] + neighborhoodSizeArray[i], 0, dataSize[i]-1 );

				lowPad[i] = (int)lowF[i] - paddedLow[i];
				highPad[i] = paddedHigh[i] - (int)highF[i];

				//lowF[i] = lowF[i] + lowPad[i];
				//highF[i] = highF[i] - highPad[i];
				//low[i] = low[i] + lowPad[i];
				//high[i] = high[i] - highPad[i];


			}// end for
			scalarField = new ITL_field_regular<SCALAR>( data[k],
														 nDim, lowF, highF,
														 lowPad, highPad, neighborhoodSizeArray );
			fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																	highF[0], highF[1], highF[2] );
			fprintf( stderr, "Padded Block boundary: %d %d %d %d %d %d\n", paddedLow[0], paddedLow[1], paddedLow[2],
																		   paddedHigh[0], paddedHigh[1], paddedHigh[2] );
			//#ifdef DEBUG_MODE
			SCALAR m = ITL_util<SCALAR>::Min( data[k], scalarField->getSize() );
			SCALAR M = ITL_util<SCALAR>::Max( data[k], scalarField->getSize() );
			printf( "Block value range: %g %g\n", m, M );
			//#endif

			if( lowPad[0] == 0 )
			{
				int enhancedSize[3];
				for( int i=0; i<3; i++ )
					enhancedSize[i] =  ( paddedHigh[i] - paddedLow[i] + 1 );
				int nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];
				FILE* dataFile = fopen( "/home/abon/block0_0.bin", "wb" );
				assert( dataFile != NULL );
				fwrite( enhancedSize, sizeof(int), 3, dataFile );
				fwrite( data[k], sizeof(float), nVertEnhanced, dataFile );
				fclose( dataFile );
			}
			if( lowPad[0] == 6 )
			{
				int enhancedSize[3];
				for( int i=0; i<3; i++ )
					enhancedSize[i] =  ( paddedHigh[i] - paddedLow[i] + 1 );
				int nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];
				FILE* dataFile = fopen( "/home/abon/block1_0.bin", "wb" );
				assert( dataFile != NULL );
				fwrite( enhancedSize, sizeof(int), 3, dataFile );
				fwrite( data[k], sizeof(float), nVertEnhanced, dataFile );
				fclose( dataFile );
			}
	
			// Create bin field
			//cout << "1" << endl;
			if( k == 0 && histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBin );

			globalEntropyComputer_scalar = new ITL_globalentropy<SCALAR>( binField, histogram, nBin );

				// Compute frequencies
				//cout << "0" << endl;
			globalEntropyComputer_scalar->computeHistogramFrequencies();

				// Get histogram frequencies
				//cout << "1" << endl;
			freqList[k] = new int[nBin];
			globalEntropyComputer_scalar->getHistogramFrequencies( freqList[k] );
			for( int i=0; i<nBin; i++ )
				fprintf( stderr, "%d, ", freqList[k][i] );
			fprintf( stderr, "\n" );

			// Initialize class that can compute entropy
			//cout << "2" << endl;
			localEntropyComputer_scalar = new ITL_localentropy<SCALAR>( binField, histogram, nBin );

			// Compute entropy
			//cout << "3" << endl;
			localEntropyComputer_scalar->computeLocalEntropyOfField( false );

			// Save local entropy field
			fprintf( stderr, "saving entropy field for block %d ...\n", k );
			blockSize[0] = (int)( highF[0] - lowF[0] + 1 );
			blockSize[1] = (int)( highF[1] - lowF[1] + 1 );
			blockSize[2] = (int)( highF[2] - lowF[2] + 1 );
			localEntropyList[k] = new float[blockSize[0]*blockSize[1]*blockSize[2]+6];
			localEntropyList[k][0] = diy_min_noghost[3*k];
			localEntropyList[k][1] = diy_min_noghost[3*k+1];
			localEntropyList[k][2] = diy_min_noghost[3*k+2];
			localEntropyList[k][3] = diy_max_noghost[3*k];
			localEntropyList[k][4] = diy_max_noghost[3*k+1];
			localEntropyList[k][5] = diy_max_noghost[3*k+2];
			memcpy( localEntropyList[k]+6,
					localEntropyComputer_scalar->getEntropyField()->getDataFull(),
					sizeof(SCALAR)*(blockSize[0]*blockSize[1]*blockSize[2]) );

			//int lowInt[3], highInt[3];
			//scalarField->getBounds( lowInt, highInt );
			//if( verboseMode == 1 ) printf( "Block Limits: %d, %d, %d, %d, %d, %d, local entropy computed.\n",
			//								lowInt[0], highInt[0],
			//								lowInt[1], highInt[1],
			//								lowInt[2], highInt[2] );

			// Clear up
			delete localEntropyComputer_scalar;
			delete binField;
			binField = NULL;
			delete scalarField;

		}
		if( fieldType == 1 )
		{
			// Initialize ITL scalar field with block data
			//cout << "0" << endl;
			for( int i=0; i<nDim; i++ )
			{
				neighborhoodSizeArray[i] = neighborhoodSize;
				paddedLow[i] = ITL_util<int>::clamp( (int)lowF[i] - neighborhoodSizeArray[i], 0, dataSize[i]-1 );
				paddedHigh[i] = ITL_util<int>::clamp( (int)highF[i] + neighborhoodSizeArray[i], 0, dataSize[i]-1 );

				lowPad[i] = (int)lowF[i] - paddedLow[i];
				highPad[i] = paddedHigh[i] - (int)highF[i];

				//lowF[i] = lowF[i] + lowPad[i];
				//highF[i] = highF[i] - highPad[i];
				//low[i] = low[i] + lowPad[i];
				//high[i] = high[i] - highPad[i];

			}// end for

			// Initialize ITL vector field with block data
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k],
														  nDim, lowF, highF,
														  lowPad, highPad, neighborhoodSizeArray );

			// Create bin field
			histMapper_vector->computeHistogramBinField_Vector( vectorField, &binField, nBin );

			// Initialize class that can compute entropy
			localEntropyComputer_vector = new ITL_localentropy<VECTOR3>( binField, histogram, nBin );

			// Compute entropy
			localEntropyComputer_vector->computeLocalEntropyOfField( false );

			// Save local entropy field
			fprintf( stderr, "saving entropy field for block %d ...\n", k );
			blockSize[0] = (int)( highF[0] - lowF[0] + 1 );
			blockSize[1] = (int)( highF[1] - lowF[1] + 1 );
			blockSize[2] = (int)( highF[2] - lowF[2] + 1 );
			localEntropyList[k] = new float[blockSize[0]*blockSize[1]*blockSize[2]+6];
			localEntropyList[k][0] = diy_min_noghost[3*k];
			localEntropyList[k][1] = diy_min_noghost[3*k+1];
			localEntropyList[k][2] = diy_min_noghost[3*k+2];
			localEntropyList[k][3] = diy_max_noghost[3*k];
			localEntropyList[k][4] = diy_max_noghost[3*k+1];
			localEntropyList[k][5] = diy_max_noghost[3*k+2];
			memcpy( localEntropyList[k]+6,
					localEntropyComputer_scalar->getEntropyField()->getDataFull(),
					sizeof(SCALAR)*(blockSize[0]*blockSize[1]*blockSize[2]) );

			//int lowInt[3], highInt[3];
			//vectorField->getBounds( lowInt, highInt );
			//if( verboseMode == 1 ) printf( "Block Limits: %d, %d, %d, %d, %d, %d, local entropy computed.\n",
			//							   lowInt[0], highInt[0],
			//							   lowInt[1], highInt[1],
			//							   lowInt[2], highInt[2] );

			// Clear up
			delete localEntropyComputer_vector;
			delete binField;
			binField = NULL;
			delete vectorField; 
		}

	}// End for loop
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Write local entropy field 
	if( verboseMode == 1 ) printf( "Writing local entropy field ...\n" );
	starttime = ITL_util<float>::startTimer();

	DIY_Write_open_all( outFile, 0 );
	DIY_Write_blocks_all( (void **)&localEntropyList[0], nblocks, NULL, 0, &CreateWriteType );
	DIY_Write_close_all();
	execTime[2] = ITL_util<float>::endTimer( starttime );


	execTime[2] = ITL_util<float>::endTimer( starttime );
	
	if( verboseMode == 1 ) printf( "%d: Read/Compute/Write Time: %f %f %f seconds\n", rank, execTime[0], execTime[1], execTime[2] );
	else printf( "%d, %f, %f, %f\n", rank, execTime[0], execTime[1], execTime[2] );


	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	delete [] diy_min_noghost;
	delete [] diy_max_noghost;
	delete [] diy_size_noghost;

	// Finalize MPI
	if( fieldType == 1 ) MPI_Type_free( &complex );
	MPI_Finalize();

}// End main




