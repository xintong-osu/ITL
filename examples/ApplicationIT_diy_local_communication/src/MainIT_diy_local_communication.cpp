/**
 * @file MainIT_diy_local.cpp
 * Application program for local entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon 
 */

#include <mpi.h>

#include "diy.h"
#include "util.hpp"
#include "assignment.hpp"
#include "blocking.hpp"
#include "io.hpp"
#include "merge.hpp"

//#include "bil.h"

#include "ITL_header.h"
#include "ITL_base.h"
#include "ITL_util.h"
#include "ITL_ioutil.h"
#include "ITL_vectormatrix.h"
#include "ITL_histogram.h"
#include "ITL_histogramconstants.h"
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

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR> *histMapper_scalar = NULL;

ITL_field_regular<int> *binField = NULL;

ITL_field_regular<SCALAR> *enhancedScalarField = NULL;
//ITL_field_regular<VECTOR3> *vectorField = NULL;
ITL_field_regular<SCALAR>** scalarFieldBlockArray = NULL;

ITL_globalentropy<SCALAR> *globalEntropyComputer_scalar = NULL;
ITL_localentropy<SCALAR> *localEntropyComputer_scalar = NULL;
//ITL_localentropy<VECTOR3> *localEntropyComputer_vector = NULL;

int dataSize[3];
int blockSize[3];
int tot_blocks = 512;

double execTime[3];
clock_t starttime, endtime;

int verboseMode = 1;

MPI_Datatype* SendItemType(int *cts, char** pts);
MPI_Datatype* RecvItemType(int *cts);
void Compute( int lid, int neighborhoodSize );
void Compute_oneblock( int lid,
					   SCALAR ***dataToSend, int *nDataPointsToSend,
					   float **remotePoint, int neighborhoodSize );
void expandBlockData( float **blockData,
					  float *lowF, float *highF,
		  	  	  	  int* lowPad, int* highPad,
		  	  	  	  int* paddedLow, int* paddedHigh,
		  	  	  	  float **receivedItems,
		  	  	  	  int nReceivedItems, float* receivedItemId,
		  	  	  	  float **enhancedData,
		  	  	  	  int neighborhoodSize );
int getRegionType( int x, int y, int z,
				   int *enhancedSize,
				   int* lowPad, int* highPad);
int getRegionIndex( int regionType,
					int nReturnedRegion, float* receivedRegionID );

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

	int dataSize[3];
	float lowF[3];
	float highF[3];
	int given[3] = {0, 0, 0};
	int neighborhoodSize = 6;
	int lowPad[3] = {0, 0, 0};
	int highPad[3] = {0, 0, 0};
	int paddedLow[3] = {0, 0, 0};
	int paddedHigh[3] = {0, 0, 0};
	int neighborhoodSizeArray[3] = {0, 0, 0};

	int fieldType = -1;
	int method = 0;


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
	//given[0] = atoi( ITL_util<float>::getArgWithName( "givenX", &argNames, &argValues ) );
	//given[1] = atoi( ITL_util<float>::getArgWithName( "givenY", &argNames, &argValues ) );
	//given[2] = atoi( ITL_util<float>::getArgWithName( "givenZ", &argNames, &argValues ) );
	tot_blocks = atoi( ITL_util<float>::getArgWithName( "nBlock", &argNames, &argValues ) );
	histogramLowEnd = atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramHighEnd = atof( ITL_util<float>::getArgWithName( "histHigh", &argNames, &argValues ) );
	verboseMode = atoi( ITL_util<float>::getArgWithName( "verbose", &argNames, &argValues ) );

	// Initialize ITL
	ITL_base::ITL_init();

	// Initialize histogram
	histogram = new ITL_histogram( patchFile, nBin );

	// Initialize data-to-histogram converter class
	histMapper_scalar = new ITL_histogrammapper<SCALAR>( histogram );

	// Allocate memory for pointers that will hold block data
	//MPI_Datatype complex;
	//if( fieldType == 1 )
	//{
	 //	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  	//	MPI_Type_commit( &complex );
	//}

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Decompose domain (with ghost layers)
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	DIY_Decompose( 0, 0, 0, given );

	// Allocate memory for pointers that will hold block data
	float* data[nblocks];
	VECTOR3* vectordata[nblocks];
	if( fieldType == 0 ) memset( data, 0, sizeof(float*) * nblocks );
	//if( fieldType == 1 ) memset( vectordata, 0, sizeof(VECTOR3*) * nblocks );

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
		//if( fieldType == 1 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(vectordata[i]));

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

	// Compute ghost layers for all blocks
	scalarFieldBlockArray = new ITL_field_regular<SCALAR>*[nblocks];
	fprintf( stderr, "Creating local fields ...\n" );
	for( int k=0; k<nblocks; k++ )
	{
		// Compute block extent (without ghost layers)
		lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
		lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
		lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

		float m = ITL_util<float>::Min( data[k], (highF[0] - lowF[0] + 1)*(highF[1] - lowF[1] + 1)*(highF[2] - lowF[2] + 1) );
		float M = ITL_util<float>::Max( data[k], (highF[0] - lowF[0] + 1)*(highF[1] - lowF[1] + 1)*(highF[2] - lowF[2] + 1) );
		fprintf( stderr, "here m M: %g %g\n", m, M );

		// Create a temporary block from original data (without ghost layers)
		scalarFieldBlockArray[k] = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );
	}
	fprintf( stderr, "Creating local fields done ...\n" );

	for( int k=0; k<nblocks; k++ )
	{
		// Create chunks of data to be sent for this block
		fprintf( stderr, "Entering compute ...\n" );
		Compute( k, neighborhoodSize );
	}

	// Scan through the blocks, compute local entropy fields
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<SCALAR> *scalarField = NULL;

	// Allocate memory for storing received ghost layers and the local entropy fields for each block
	float* localEntropyList[nblocks];
	float* enhancedData[nblocks];
	void ***receivedItems = new void**[nblocks]; // received items from neighbors (generic pointer to a byte, not a string)
	int *num_items = new int[nblocks]; // number of received items in each block
	float *returnedItemId = new float[26]; // number of received items in each block
  	int *freqList[nblocks];

	// exchange neighbors
  	fprintf( stderr, "Starting Neighbor exchange ...\n" );
	DIY_Exchange_neighbors( receivedItems, num_items, 1.0, &RecvItemType, &SendItemType );
	fprintf( stderr, "Neighbor exchange done ...\n" );

	for( int k=0; k<nblocks; k++ )
	{
		if( fieldType == 0 )
		{
			// Compute block extent (without ghost layers)
			lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
			lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
			lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];
			for( int i=0; i<nDim; i++ )
			{
				neighborhoodSizeArray[i] = neighborhoodSize;
				paddedLow[i] = ITL_util<int>::clamp( (int)lowF[i] - neighborhoodSizeArray[i], 0, dataSize[i]-1 );
				paddedHigh[i] = ITL_util<int>::clamp( (int)highF[i] + neighborhoodSizeArray[i], 0, dataSize[i]-1 );

				lowPad[i] = (int)lowF[i] - paddedLow[i];
				highPad[i] = paddedHigh[i] - (int)highF[i];

			}// end for
			fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																	highF[0], highF[1], highF[2] );
			fprintf( stderr, "Padded Block boundary: %d %d %d %d %d %d\n", paddedLow[0], paddedLow[1], paddedLow[2],
																		   paddedHigh[0], paddedHigh[1], paddedHigh[2] );



			// Find order of the returned items
			fprintf( stderr, "Number of items received by block %d: %d \n", k, num_items[k] );
			for( int i=0; i<num_items[k]; i++ )
			{
				assert( receivedItems[k][i] != NULL );
				returnedItemId[i] = ((float*)receivedItems[k][i])[0];
				fprintf( stderr, "%d: %g\n", i, returnedItemId[i] );
			}

			// Expand block data with ghost layers
			fprintf( stderr, "expanding block data for block %d ...\n", k );
			expandBlockData( &data[k],
							 lowF, highF,
							 lowPad, highPad,
							 paddedLow,paddedHigh,
							 (float**)receivedItems[k],
							 num_items[k], returnedItemId,
							 &enhancedData[k], neighborhoodSize );

			// Initialize ITL scalar field with block data
			fprintf( stderr, "creating scalar field for block %d ...\n", k );
			scalarField = new ITL_field_regular<SCALAR>( enhancedData[k], nDim, lowF, highF,
														 lowPad, highPad, neighborhoodSizeArray );

			// Create bin field
			fprintf( stderr, "creating mapper for block %d ...\n", k );
			if( k == 0 && histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBin );

			int enhancedBlockSize = (paddedHigh[0] - paddedLow[0] + 1)*(paddedHigh[1] - paddedLow[1] + 1)*(paddedHigh[2] - paddedLow[2] + 1);
			int binm = ITL_util<int>::Min( binField->getDataFull(), enhancedBlockSize  );
			int binM = ITL_util<int>::Max( binField->getDataFull(), enhancedBlockSize );
			fprintf( stderr, "binField m M: %d %d\n", binm, binM );

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
			fprintf( stderr, "creating entropy class for block %d ...\n", k );
			localEntropyComputer_scalar = new ITL_localentropy<SCALAR>( binField, histogram, nBin );


			// Compute entropy
			fprintf( stderr, "computing entropy for block %d ...\n", k );
			localEntropyComputer_scalar->computeLocalEntropyOfField( false );

			// Save local entropy field
			fprintf( stderr, "saving entropy field for block %d ...\n", k );
			blockSize[0] = (int)( highF[0] - lowF[0] + 1 );
			blockSize[1] = (int)( highF[1] - lowF[1] + 1 );
			blockSize[2] = (int)( highF[2] - lowF[2] + 1 );
			localEntropyList[k] = new float[blockSize[0]*blockSize[1]*blockSize[2]+6];
			localEntropyList[k][0] = diy_min[3*k];
			localEntropyList[k][1] = diy_min[3*k+1];
			localEntropyList[k][2] = diy_min[3*k+2];
			localEntropyList[k][3] = diy_max[3*k];
			localEntropyList[k][4] = diy_max[3*k+1];
			localEntropyList[k][5] = diy_max[3*k+2];
			memcpy( localEntropyList[k]+6,
					localEntropyComputer_scalar->getEntropyField()->getDataFull(),
					sizeof(SCALAR)*(blockSize[0]*blockSize[1]*blockSize[2]) );

			SCALAR em = ITL_util<float>::Min( localEntropyComputer_scalar->getEntropyField()->getDataFull(), blockSize[0]*blockSize[1]*blockSize[2]  );
			SCALAR eM = ITL_util<float>::Max( localEntropyComputer_scalar->getEntropyField()->getDataFull(), blockSize[0]*blockSize[1]*blockSize[2] );
			fprintf( stderr, "enField m M: %g %g\n", em, eM );

			// Clear up
			delete localEntropyComputer_scalar;
			delete binField;
			binField = NULL;
			delete scalarField;
		}
		/*
		if( fieldType == 1 )
		{

			// Initialize ITL vector field with block data
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

			// Initialize class that can compute entropy
			localEntropyComputer_vector = new ITL_localentropy<VECTOR3>( vectorField );

			// Create bin field
			localEntropyComputer_vector->computeHistogramBinField( "vector", nBin );

			// Compute entropy
			localEntropyComputer_vector->computeEntropyOfField( nBin, false );

			// Save entropy
			localEntropyList[k] = new float[vectorField->grid->nVertices];
			memcpy( localEntropyList[k], localEntropyComputer_vector->getEntropyField(), sizeof(SCALAR)*vectorField->grid->nVertices );
			if( verboseMode == 1 ) printf( "Block Limits: %d, %d, %d, %d, %d, %d, local entropy computed.\n", vectorField->grid->lowInt[0],
															vectorField->grid->highInt[0],
															vectorField->grid->lowInt[1],
															vectorField->grid->highInt[1],
															vectorField->grid->lowInt[2],
															vectorField->grid->highInt[2] );

			// Clear up
			delete localEntropyComputer_vector;
			delete vectorField;
		}
		*/

	}// End for loop

	// flush any remaining messages
	// if multiple rounds of compute / exchange neighbors, call FlushNeighbors
	// once for each time block, after the rounds complete
	DIY_Flush_neighbors( receivedItems, num_items, &RecvItemType );

	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Write local entropy field
	if( verboseMode == 1 ) printf( "Writing local entropy field ...\n" );
	starttime = ITL_util<float>::startTimer();

	DIY_Write_open_all( outFile, 0 );
	DIY_Write_blocks_all( (void **)&localEntropyList[0], nblocks, NULL, 0, &CreateWriteType );
	DIY_Write_close_all();
	execTime[2] = ITL_util<float>::endTimer( starttime );

	if( verboseMode == 1 ) printf( "%d: Read/Compute/Write Time: %f %f %f seconds\n", rank, execTime[0], execTime[1], execTime[2] );
	else printf( "%d, %f, %f, %f\n", rank, execTime[0], execTime[1], execTime[2] );

	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;

	// Finalize MPI
	//if( fieldType == 1 ) MPI_Type_free( &complex );
	MPI_Finalize();

}// End main




//
// computation for my local block number lid (local ID)
// Here I assume x spans from left to right, y spans from front to back and z spans from bottom to top
//
void Compute( int lid, int neighborhoodSize )
{
	int maxDimLength = 64;
	int maxGhostSize = neighborhoodSize * maxDimLength * maxDimLength;
	int lowBoundary[3];
	int highBoundary[3];
	int nItemToSendPerBlock = 26;
	int nDataPointsToSend[nItemToSendPerBlock];
	SCALAR** dataToSend = new SCALAR*[nItemToSendPerBlock];
	//float** remotePoint = new float*[nItemToSendPerBlock];
	for( int i=0; i<nItemToSendPerBlock; i++ )
	{
		dataToSend[i] = NULL;
		//remotePoint[i] = new float[3];
	}
	unsigned char dirs[26] = { DIY_X0,
							 DIY_Y0,
							 DIY_X1,
							 DIY_Y1,
							 DIY_Z0,
							 DIY_Z1,
							 DIY_X0 | DIY_Z0,
							 DIY_Y0 | DIY_Z0,
							 DIY_X1 | DIY_Z0,
							 DIY_Y1 | DIY_Z0,
							 DIY_X0 | DIY_Z1,
							 DIY_Y0 | DIY_Z1,
							 DIY_X1 | DIY_Z1,
							 DIY_Y1 | DIY_Z1,
							 DIY_X0 | DIY_Y0,
							 DIY_X1 | DIY_Y0,
							 DIY_X1 | DIY_Y1,
							 DIY_X0 | DIY_Y1,
							 DIY_X0 | DIY_Y0 | DIY_Z0,
							 DIY_X1 | DIY_Y0 | DIY_Z0,
							 DIY_X1 | DIY_Y1 | DIY_Z0,
							 DIY_X0 | DIY_Y1 | DIY_Z0,
							 DIY_X0 | DIY_Y0 | DIY_Z1,
							 DIY_X1 | DIY_Y0 | DIY_Z1,
							 DIY_X1 | DIY_Y1 | DIY_Z1,
							 DIY_X0 | DIY_Y1 | DIY_Z1 };

	// local computation here, producing an item that needs to be sent to 
  	// a neighboring block
	fprintf( stderr, "Entering compute one block for block %d ..\n", lid );

	// 6 ghost faces
	// First we iterate through the side faces in a counterclockwise manner
	// Left (id: 1)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
    //fprintf( stderr, "got bounds %d\n", lid );
	highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
	nDataPointsToSend[0] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[0] );
	}
	assert( nDataPointsToSend[0] > 0 );
    //dataToSend[0] = new SCALAR[1+nDataPointsToSend[0]];
    dataToSend[0] = new SCALAR[1+maxGhostSize];
    dataToSend[0][0] = 1;
    //fprintf( stderr, "allocated %d\n", lid );
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[0][1] );
	float m = ITL_util<float>::Min( &dataToSend[0][1], nDataPointsToSend[0] );
	float M = ITL_util<float>::Max( &dataToSend[0][1], nDataPointsToSend[0] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[0], m, M );
	//fprintf( stderr, "got data %d\n", lid );
	//remotePoint[0][0] = lowBoundary[0] - 1;
	//remotePoint[0][1] = lowBoundary[1] + 1;
	//remotePoint[0][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 0 complete\n" );

	// Front (id: 2)
	scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	nDataPointsToSend[1] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[1] );
	}
	//dataToSend[1] = new SCALAR[1+nDataPointsToSend[1]];
	dataToSend[1] = new SCALAR[1+maxGhostSize];
	dataToSend[1][0] = 2;
	//fprintf( stderr, "allocated %d\n", lid );
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[1][1] );
	m = ITL_util<float>::Min( &dataToSend[1][1], nDataPointsToSend[1] );
	M = ITL_util<float>::Max( &dataToSend[1][1], nDataPointsToSend[1] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[1], m, M );
	//fprintf( stderr, "got data %d\n", lid );
	//remotePoint[1][0] = lowBoundary[0] + 1;
	//remotePoint[1][1] = lowBoundary[1] + 1;
	//remotePoint[1][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 1 complete\n" );

	// Right (id: 3)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
 	nDataPointsToSend[2] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[2] );
	}
    //dataToSend[2] = new SCALAR[1+nDataPointsToSend[2]];
    dataToSend[2] = new SCALAR[1+maxGhostSize];
    dataToSend[2][0] = 3;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[2][1] );
	m = ITL_util<float>::Min( &dataToSend[2][1], nDataPointsToSend[2] );
	M = ITL_util<float>::Max( &dataToSend[2][1], nDataPointsToSend[2] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[2], m, M );

	//remotePoint[2][0] = highBoundary[0] + 1;
	//remotePoint[2][1] = highBoundary[1] - 1;
	//remotePoint[2][2] = highBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 2 complete\n" );

	// Back (id: 4)
	scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	nDataPointsToSend[3] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[3] );
	}
    //dataToSend[3] = new SCALAR[1+nDataPointsToSend[3]];
    dataToSend[3] = new SCALAR[1+maxGhostSize];
    dataToSend[3][0] = 4;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[3][1] );
	m = ITL_util<float>::Min( &dataToSend[3][1], nDataPointsToSend[3] );
	M = ITL_util<float>::Max( &dataToSend[3][1], nDataPointsToSend[3] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[3], m, M );

	//remotePoint[3][0] = highBoundary[0] - 1;
	//remotePoint[3][1] = highBoundary[1] - 1;
	//remotePoint[3][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 3 complete\n" );

	// Now we do the bottom and the top face
	// Bottom (id: 5)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[4] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[4] );
	}
    //dataToSend[4] = new SCALAR[1+nDataPointsToSend[4]];
    dataToSend[4] = new SCALAR[1+maxGhostSize];
    dataToSend[4][0] = 5;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[4][1] );
	m = ITL_util<float>::Min( &dataToSend[4][1], nDataPointsToSend[4] );
	M = ITL_util<float>::Max( &dataToSend[4][1], nDataPointsToSend[4] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[4], m, M );

	//remotePoint[4][0] = highBoundary[0] - 1;
	//remotePoint[4][1] = highBoundary[1] + 1;
	//remotePoint[4][2] = highBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 4 complete\n" );

	// Top (id: 6)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[5] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[5] );
	}
    //dataToSend[5] = new SCALAR[1+nDataPointsToSend[5]];
    dataToSend[5] = new SCALAR[1+maxGhostSize];
    dataToSend[5][0] = 6;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[5][1] );
	m = ITL_util<float>::Min( &dataToSend[5][1], nDataPointsToSend[5] );
	M = ITL_util<float>::Max( &dataToSend[5][1], nDataPointsToSend[5] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[5], m, M );

	//remotePoint[5][0] = lowBoundary[0] + 1;
	//remotePoint[5][1] = lowBoundary[1] - 1;
	//remotePoint[5][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 5 complete\n" );
	
	// 12 ghost bars grazing the edges 
	// First we do the bottom 4 edges in a counterclockwise manner	
	// Bottom left (id: 7)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[6] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[6] );
	}
    //dataToSend[6] = new SCALAR[1+nDataPointsToSend[6]];
    dataToSend[6] = new SCALAR[1+maxGhostSize];
    dataToSend[6][0] = 7;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[6][1] );
	m = ITL_util<float>::Min( &dataToSend[6][1], nDataPointsToSend[6] );
	M = ITL_util<float>::Max( &dataToSend[6][1], nDataPointsToSend[6] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[6], m, M );
	//remotePoint[6][0] = lowBoundary[0] - 1;
	//remotePoint[6][1] = lowBoundary[1] - 1;
	//remotePoint[6][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 6 complete\n" );

	// Bottom front (id: 8)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[7] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[7] );
	}
    //dataToSend[7] = new SCALAR[1+nDataPointsToSend[7]];
    dataToSend[7] = new SCALAR[1+maxGhostSize];
    dataToSend[7][0] = 8;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[7][1] );
	m = ITL_util<float>::Min( &dataToSend[7][1], nDataPointsToSend[7] );
	M = ITL_util<float>::Max( &dataToSend[7][1], nDataPointsToSend[7] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[7], m, M );
	//remotePoint[7][0] = lowBoundary[0] + 1;
	//remotePoint[7][1] = lowBoundary[1] - 1;
	//remotePoint[7][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 7 complete\n" );

	// Bottom right (id: 9)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
   	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[8] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[8] );
	}
    //dataToSend[8] = new SCALAR[1+nDataPointsToSend[8]];
    dataToSend[8] = new SCALAR[1+maxGhostSize];
    dataToSend[8][0] = 9;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[8][1] );
	m = ITL_util<float>::Min( &dataToSend[8][1], nDataPointsToSend[8] );
	M = ITL_util<float>::Max( &dataToSend[8][1], nDataPointsToSend[8] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[8], m, M );
	//remotePoint[8][0] = highBoundary[0] + 1;
	//remotePoint[8][1] = lowBoundary[1] - 1;
	//remotePoint[8][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 8 complete\n" );

	// Bottom back (id: 10)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
  	highBoundary[2] = lowBoundary[2] + neighborhoodSize + 1;
	nDataPointsToSend[9] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[9] );
	}
    //dataToSend[9] = new SCALAR[1+nDataPointsToSend[9]];
    dataToSend[9] = new SCALAR[1+maxGhostSize];
    dataToSend[9][0] = 10;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[9][1] );
	m = ITL_util<float>::Min( &dataToSend[9][1], nDataPointsToSend[9] );
	M = ITL_util<float>::Max( &dataToSend[9][1], nDataPointsToSend[9] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[9], m, M );
	//remotePoint[9][0] = lowBoundary[0] + 1;
	//remotePoint[9][1] = lowBoundary[1] - 1;
	//remotePoint[9][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 9 complete\n" );

	// Next we do the top 4 edges in a counterclockwise manner	
	// Top left (id: 11)	
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
  	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	highBoundary[0] = lowBoundary[0] + neighborhoodSize + 1;
	nDataPointsToSend[10] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[10] );
	}
    //dataToSend[10] = new SCALAR[1+nDataPointsToSend[10]];
    dataToSend[10] = new SCALAR[1+maxGhostSize];
    dataToSend[10][0] = 11;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[10][1] );
	m = ITL_util<float>::Min( &dataToSend[10][1], nDataPointsToSend[10] );
	M = ITL_util<float>::Max( &dataToSend[10][1], nDataPointsToSend[10] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[10], m, M );
	//remotePoint[10][0] = lowBoundary[0] - 1;
	//remotePoint[10][1] = highBoundary[1] + 1;
	//remotePoint[10][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 10 complete\n" );

	// Top front (id: 12)	
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
  	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	nDataPointsToSend[11] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[11] );
	}
    //dataToSend[11] = new SCALAR[1+nDataPointsToSend[11]];
    dataToSend[11] = new SCALAR[1+maxGhostSize];
    dataToSend[11][0] = 12;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[11][1] );
	m = ITL_util<float>::Min( &dataToSend[11][1], nDataPointsToSend[11] );
	M = ITL_util<float>::Max( &dataToSend[11][1], nDataPointsToSend[11] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[11], m, M );
	//remotePoint[11][0] = lowBoundary[0] + 1;
	//remotePoint[11][1] = highBoundary[1] + 1;
	//remotePoint[11][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 11 complete\n" );

	// Top right (id: 13)	
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
  	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[12] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[12] );
	}
    //dataToSend[12] = new SCALAR[1+nDataPointsToSend[12]];
    dataToSend[12] = new SCALAR[1+maxGhostSize];
    dataToSend[12][0] = 13;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[12][1] );
	m = ITL_util<float>::Min( &dataToSend[12][1], nDataPointsToSend[12] );
	M = ITL_util<float>::Max( &dataToSend[12][1], nDataPointsToSend[12] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[12], m, M );

	//remotePoint[12][0] = highBoundary[0] + 1;
	//remotePoint[12][1] = highBoundary[1] + 1;
	//remotePoint[12][2] = lowBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 12 complete\n" );

	// Top back (id: 14)	
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[13] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[13] );
	}
    //dataToSend[13] = new SCALAR[1+nDataPointsToSend[13]];
    dataToSend[13] = new SCALAR[1+maxGhostSize];
    dataToSend[13][0] = 14;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[13][1] );
	m = ITL_util<float>::Min( &dataToSend[13][1], nDataPointsToSend[13] );
	M = ITL_util<float>::Max( &dataToSend[13][1], nDataPointsToSend[13] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[13], m, M );
	//remotePoint[13][0] = lowBoundary[0] + 1;
	//remotePoint[13][1] = highBoundary[1] + 1;
	//remotePoint[13][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 13 complete\n" );

	// Next we do the standing 4 edges in a counterclockwise manner	
	// Standing front left (id: 15) 
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[0] = lowBoundary[0] + neighborhoodSize -1;
  	highBoundary[1] = lowBoundary[1] + neighborhoodSize -1;
	nDataPointsToSend[14] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[14] );
	}
    //dataToSend[14] = new SCALAR[1+nDataPointsToSend[14]];
    dataToSend[14] = new SCALAR[1+maxGhostSize];
    dataToSend[14][0] = 15;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[14][1] );
	m = ITL_util<float>::Min( &dataToSend[14][1], nDataPointsToSend[14] );
	M = ITL_util<float>::Max( &dataToSend[14][1], nDataPointsToSend[14] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[14], m, M );
	//remotePoint[14][0] = lowBoundary[0] - 1;
	//remotePoint[14][1] = lowBoundary[1] + 1;
	//remotePoint[14][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 14 complete\n" );

	// Standing front right (id: 16)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
 	highBoundary[1] = lowBoundary[1] + neighborhoodSize -1;
	nDataPointsToSend[15] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[15] );
	}
    //dataToSend[15] = new SCALAR[1+nDataPointsToSend[15]];
    dataToSend[15] = new SCALAR[1+maxGhostSize];
    dataToSend[15][0] = 16;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[15][1] );
	m = ITL_util<float>::Min( &dataToSend[15][1], nDataPointsToSend[15] );
	M = ITL_util<float>::Max( &dataToSend[15][1], nDataPointsToSend[15] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[15], m, M );

	//remotePoint[15][0] = highBoundary[0] + 1;
	//remotePoint[15][1] = lowBoundary[1] + 1;
	//remotePoint[15][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 15 complete\n" );

	// Standing back right (id: 17)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
 	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	nDataPointsToSend[16] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[16] );
	}
    //dataToSend[16] = new SCALAR[1+nDataPointsToSend[16]];
    dataToSend[16] = new SCALAR[1+maxGhostSize];
    dataToSend[16][0] = 17;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[16][1] );
	m = ITL_util<float>::Min( &dataToSend[16][1], nDataPointsToSend[16] );
	M = ITL_util<float>::Max( &dataToSend[16][1], nDataPointsToSend[16] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[16], m, M );
	//remotePoint[16][0] = highBoundary[0] + 1;
	//remotePoint[16][1] = lowBoundary[1] + 1;
	//remotePoint[16][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 16 complete\n" );

	// Standing back left (id: 18)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
    lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	nDataPointsToSend[17] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[17] );
	}
    //dataToSend[17] = new SCALAR[1+nDataPointsToSend[17]];
    dataToSend[17] = new SCALAR[1+maxGhostSize];
    dataToSend[17][0] = 18;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[17][1] );
	m = ITL_util<float>::Min( &dataToSend[17][1], nDataPointsToSend[17] );
	M = ITL_util<float>::Max( &dataToSend[17][1], nDataPointsToSend[17] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[17], m, M );
	//remotePoint[17][0] = lowBoundary[0] - 1;
	//remotePoint[17][1] = lowBoundary[1] + 1;
	//remotePoint[17][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 17 complete\n" );

	// 8 corners
	// First we do the bottom 4 corners in a counterclockwise manner	
	// Bottom front left (id: 19)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[18] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[18] );
	}
    //dataToSend[18] = new SCALAR[1+nDataPointsToSend[18]];
    dataToSend[18] = new SCALAR[1+maxGhostSize];
    dataToSend[18][0] = 19;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[18][1] );
	m = ITL_util<float>::Min( &dataToSend[18][1], nDataPointsToSend[18] );
	M = ITL_util<float>::Max( &dataToSend[18][1], nDataPointsToSend[18] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[18], m, M );
	//remotePoint[18][0] = lowBoundary[0] - 1;
	//remotePoint[18][1] = lowBoundary[1] - 1;
	//remotePoint[18][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 18 complete\n" );

	// Bottom front right (id: 20)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[19] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[19] );
	}
    //dataToSend[19] = new SCALAR[1+nDataPointsToSend[19]];
    dataToSend[19] = new SCALAR[1+maxGhostSize];
    dataToSend[19][0] = 20;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[19][1] );
	m = ITL_util<float>::Min( &dataToSend[19][1], nDataPointsToSend[19] );
	M = ITL_util<float>::Max( &dataToSend[19][1], nDataPointsToSend[19] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[19], m, M );
	//remotePoint[19][0] = highBoundary[0] - 1;
	//remotePoint[19][1] = highBoundary[1] + 1;
	//remotePoint[19][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 19 complete\n" );

	// Bottom back right (id: 21)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
 	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	//highBoundary[2] = highBoundary[2];
	nDataPointsToSend[20] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[20] );
	}
    //dataToSend[20] = new SCALAR[1+nDataPointsToSend[20]];
    dataToSend[20] = new SCALAR[1+maxGhostSize];
    dataToSend[20][0] = 21;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[20][1] );
	m = ITL_util<float>::Min( &dataToSend[20][1], nDataPointsToSend[20] );
	M = ITL_util<float>::Max( &dataToSend[20][1], nDataPointsToSend[20] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[20], m, M );
	//remotePoint[20][0] = highBoundary[0] + 1;
	//remotePoint[20][1] = lowBoundary[1] - 1;
	//remotePoint[20][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 20 complete\n" );

	// Bottom back left (id: 22)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
	nDataPointsToSend[21] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[21] );
	}
    //dataToSend[21] = new SCALAR[1+nDataPointsToSend[21]];
    dataToSend[21] = new SCALAR[1+maxGhostSize];
    dataToSend[21][0] = 22;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[21][1] );
	m = ITL_util<float>::Min( &dataToSend[21][1], nDataPointsToSend[21] );
	M = ITL_util<float>::Max( &dataToSend[21][1], nDataPointsToSend[21] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[21], m, M );
	//remotePoint[21][0] = lowBoundary[0] - 1;
	//remotePoint[21][1] = lowBoundary[1] - 1;
	//remotePoint[21][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 21 complete\n" );

	// Next we do the top 4 corners in a counterclockwise manner	
	// Top front left (id: 23)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
    highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
    lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[22] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[22] );
	}
    //dataToSend[22] = new SCALAR[1+nDataPointsToSend[22]];
    dataToSend[22] = new SCALAR[1+maxGhostSize];
    dataToSend[22][0] = 23;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[22][1] );
	m = ITL_util<float>::Min( &dataToSend[22][1], nDataPointsToSend[22] );
	M = ITL_util<float>::Max( &dataToSend[22][1], nDataPointsToSend[22] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[22], m, M );
	//remotePoint[22][0] = lowBoundary[0] - 1;
	//remotePoint[22][1] = highBoundary[1] + 1;
	//remotePoint[22][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 22 complete\n" );

	// Top front right (id: 24)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[23] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[23] );
	}
    //dataToSend[23] = new SCALAR[1+nDataPointsToSend[23]];
    dataToSend[23] = new SCALAR[1+maxGhostSize];
    dataToSend[23][0] = 24;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[23][1] );
	m = ITL_util<float>::Min( &dataToSend[23][1], nDataPointsToSend[23] );
	M = ITL_util<float>::Max( &dataToSend[23][1], nDataPointsToSend[23] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[23], m, M );
	//remotePoint[23][0] = highBoundary[0] + 1;
	//remotePoint[23][1] = highBoundary[1] + 1;
	//remotePoint[23][2] = lowBoundary[2] - 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 23 complete\n" );

	// Top back right (id: 25)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
	lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[24] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[24] );
	}
    //dataToSend[24] = new SCALAR[1+nDataPointsToSend[24]];
    dataToSend[24] = new SCALAR[1+maxGhostSize];
    dataToSend[24][0] = 25;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[24][1] );
	m = ITL_util<float>::Min( &dataToSend[24][1], nDataPointsToSend[24] );
	M = ITL_util<float>::Max( &dataToSend[24][1], nDataPointsToSend[24] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[24], m, M );
	//remotePoint[24][0] = highBoundary[0] + 1;
	//remotePoint[24][1] = highBoundary[1] + 1;
	//remotePoint[24][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 24 complete\n" );

	// Top back left (id: 26)
    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );
    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
    lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
	nDataPointsToSend[25] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
						   ( highBoundary[1] - lowBoundary[1] + 1 ) *
						   ( highBoundary[2] - lowBoundary[2] + 1 );
	if( verboseMode == 1 )
	{
		fprintf( stderr, "boundary: %d %d %d %d %d %d\n", lowBoundary[0], lowBoundary[1], lowBoundary[2],
														  highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "pad size %d\n", nDataPointsToSend[25] );
	}
    //dataToSend[25] = new SCALAR[1+nDataPointsToSend[25]];
    dataToSend[25] = new SCALAR[1+maxGhostSize];
    dataToSend[25][0] = 26;
	scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[25][1] );
	m= ITL_util<float>::Min( &dataToSend[25][1], nDataPointsToSend[25] );
	M = ITL_util<float>::Max( &dataToSend[25][1], nDataPointsToSend[25] );
	fprintf( stderr, "packet size and range: %d %g %g\n", nDataPointsToSend[25], m, M );
	//remotePoint[25][0] = lowBoundary[0] - 1;
	//remotePoint[25][1] = highBoundary[1] + 1;
	//remotePoint[25][2] = highBoundary[2] + 1;
	if( verboseMode == 1 ) fprintf( stderr, "Item 25 complete\n" );

   	fprintf( stderr, "Out from compute one block for block %d ..\n", lid );

	// enqueue the item for sending to neighbor
	// in this example the point by which the neighbor is identified (5th arg.)
	// is the same as the item (2nd arg.), but this need not be the case
	// because the item can be any generic data

	for( int j=0; j<nItemToSendPerBlock; j++ )
	{
		//DIY_Enqueue_item_points( lid, dataToSend[j], NULL,
		//						 nDataPointsToSend[j] * sizeof(float),
		//						 remotePoint[j], 1, NULL );
		DIY_Enqueue_item_dirs( lid, &dataToSend[j][0], NULL,
							   (nDataPointsToSend[j]+1) * sizeof(float),
							   &dirs[j], 1, NULL );
	}// end inner for
	fprintf( stderr, "enqueue item done for block %d ..\n", lid );

}// end for

// Here I assume x spans from left to right, y spans from bottom to top and z spans from front to back
void expandBlockData( float **blockData, float *lowF, float *highF,
					  int* lowPad, int* highPad,
					  int* paddedLow, int* paddedHigh,
					  float **receivedItems,
					  int nReceivedItems, float* receivedItemId,
					  float **enhancedData,
					  int neighborhoodSize )
{
	int blocksize[3], enhancedSize[3];
	int recvdItemIndex[nReceivedItems];
	for( int i=0; i<nReceivedItems; i++ )
		recvdItemIndex[i] = 1;
	fprintf( stderr, "%g\n", receivedItems[0][0] );
	
	// Allocate memory for enhanced block
	for( int i=0; i<3; i++ )
	{
		blocksize[i] = (int)( highF[i] - lowF[i] + 1 );
		enhancedSize[i] = ( paddedHigh[i] - paddedLow[i] + 1 );
	}
	int nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];
	(*enhancedData) = new float[nVertEnhanced];

	//if( verboseMode == 1 )
	fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
															highF[0], highF[1], highF[2] );
	fprintf( stderr, "Pad size: %d %d %d %d %d %d\n", lowPad[0], lowPad[1], lowPad[2],
													  highPad[0], highPad[1], highPad[2] );
	fprintf( stderr, "Enhanced size: %d %d %d\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );

	float m = ITL_util<float>::Min( (*blockData), blocksize[0]*blocksize[1]*blocksize[2] );
	float M = ITL_util<float>::Max( (*blockData), blocksize[0]*blocksize[1]*blocksize[2] );
	fprintf( stderr, "m M: %g %g\n", m, M );

	m = ITL_util<float>::Min( (*receivedItems) + 1, 6*64*64 );
	M = ITL_util<float>::Max( (*receivedItems) + 1, 6*64*64 );
	fprintf( stderr, "rm rM: %g %g\n", m, M );

	int indexid = 0;
	int ownDataIndex = 0;
	int mirrorDataCount = 0;
	int mirrorx, mirrory, mirrorz, mirrorIndex;
	for( int z=0; z<enhancedSize[2]; z++ )
	{	
		for( int y=0; y<enhancedSize[1]; y++ )
		{
			for( int x=0; x<enhancedSize[0]; x++ )
			{
				//fprintf( stderr, "%d %d %d\n", x, y, z );

				if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
					y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
					z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
				// Inside its own area
				{
					//fprintf( stderr, "Data value: %g\n", (*blockData)[ownDataIndex] );
					(*enhancedData)[indexid] = (*blockData)[ownDataIndex];
					ownDataIndex++;

					//if( x > enhancedSize[0]-highPad[0] )
					//{
					//	fprintf( stderr, "%d %d %d\n", x, y, z );
					//	fprintf( stderr, "Enhanced size: %d %d %d\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );
					//	fprintf( stderr, "Pad size: %d %d %d %d %d %d\n", lowPad[0], lowPad[1], lowPad[2],
					//			 highPad[0], highPad[1], highPad[2] );
					//	exit(0);
					//}
				}
				else
				// Inside one of the ghost layers
				{
					// Get Region Type
					int regionType = getRegionType( x, y, z, enhancedSize, lowPad, highPad );

					// Get index of region
					int itemIndex = getRegionIndex( regionType, nReceivedItems, receivedItemId );
					//fprintf( stderr, "region type and index: %d %d\n", regionType, itemIndex );

					if( itemIndex == -1 )
					{
						// Mirror data
						//mirrorx = ITL_util<int>::mirror( x, lowPad[0], lowPad[0] + blocksize[0]-1 ) - lowPad[0];
						//mirrory = ITL_util<int>::mirror( y, lowPad[1], lowPad[1] + blocksize[1]-1 ) - lowPad[1];
						//mirrorz = ITL_util<int>::mirror( z, lowPad[2], lowPad[2] + blocksize[2]-1 ) - lowPad[2];
						//mirrorIndex = mirrorz * blocksize[0] * blocksize[1] +
						//			  mirrory * blocksize[0] +
						//			  mirrorx;
						//fprintf( stderr, "mirrored: %d %d %d\n", mirrorx, mirrory, mirrorz );

						//(*enhancedData)[indexid] = (*blockData)[mirrorIndex];
						mirrorDataCount ++;
					}
					else
					{
						//if( regionType != 1 && regionType != 3 )
						//{
						//	fprintf( stderr, "%d %d %d %d\n", x, y, z, regionType );
						//	fprintf( stderr, "Enhanced size: %d %d %d\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );
						//	fprintf( stderr, "Pad size: %d %d %d %d %d %d\n", lowPad[0], lowPad[1], lowPad[2],
						//						highPad[0], highPad[1], highPad[2] );
						//	exit(0);
						//}

						(*enhancedData)[indexid] = receivedItems[itemIndex][recvdItemIndex[itemIndex]];
						recvdItemIndex[itemIndex]++;
					}
				}

				indexid++;
 
			}// end for z
		}// end for y
	}// end for z

	if( lowPad[0] == 0 )
	{
		FILE* dataFile = fopen( "/home/abon/block0_1.bin", "wb" );
		assert( dataFile != NULL );
		fwrite( enhancedSize, sizeof(int), 3, dataFile );
		fwrite( (*enhancedData), sizeof(float), nVertEnhanced, dataFile );
		fclose( dataFile );
	}
	if( lowPad[0] == 6 )
	{
		FILE* dataFile = fopen( "/home/abon/block1_1.bin", "wb" );
		assert( dataFile != NULL );
		fwrite( enhancedSize, sizeof(int), 3, dataFile );
		fwrite( (*enhancedData), sizeof(float), nVertEnhanced, dataFile );
		fclose( dataFile );
	}


	
	for( int i=0; i<nReceivedItems; i++ )
		fprintf( stderr, "%d -> %d\n", i, recvdItemIndex[i] );
	fprintf( stderr, "inside points: %d\n", ownDataIndex );
	fprintf( stderr, "mirror points: %d\n", mirrorDataCount );


	m = ITL_util<float>::Min( (*enhancedData), enhancedSize[0]*enhancedSize[1]*enhancedSize[2] );
	M = ITL_util<float>::Max( (*enhancedData), enhancedSize[0]*enhancedSize[1]*enhancedSize[2] );
	fprintf( stderr, "xm xM: %g %g\n", m, M );

}// end function 

// Here I assume x spans from left to right, y spans from front to back and z spans from bottom to top
int
getRegionType( int x, int y, int z, int *enhancedSize, int* lowPad, int* highPad )
{
	// Determine if point in any of the 6 ghost faces
	// left	: corresponds to right
	if( x>=0 && x<lowPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 3;
	// front: corresponds to back
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 4;
	// right: corresponds to left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 1;
	// back: corresponds to front
	if( x>=lowPad[0] && x<=enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=lowPad[2] && z<=enhancedSize[2]-highPad[2] )
		return 2;
	// bottom: corresponds to top
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=0 && z<lowPad[2] )
		return 6;
	// top: corresponds to bottom
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 5;
	 
	// Determine if point in any of the 12 ghost edge-grazing bars
	// bottom left	: corresponds to top right
	if( x>=0 && x<lowPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=0 && z<lowPad[2] )
		return 12;
	// bottom front: corresponds to top back
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=0 && z<lowPad[2] )
		return 13;
	// bottom right: corresponds to top left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=0 && z<lowPad[2] )
		return 14;
	// bottom back: corresponds to top front
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=0 && z<lowPad[2] )
		return 11;
	/////////////////////////////////////////
	// top left	: corresponds to bottom right
	if( x>=0 && x<lowPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 8;
	// top front: corresponds to bottom back
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=0 && y<highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 9;
	// top right: corresponds to bottom left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 10;
	// top back: corresponds to bottom front
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 7;
	//////////////////////////////////////////////////////////
	// standing front left: corresponds to standing back right
	if( x>=0 && x<lowPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 17;
	// standing front right: corresponds to standing back left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=0 && y<lowPad[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 18;
	// standing back right: corresponds to front left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 15;
	// standing back left: corresponds to front front right
	if( x>=0 && x<lowPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
		return 16;

	// Determine if point in any of the 8 ghost corners		
	// bottom front left: corresponds to top back right
	if( x>=0 && x<lowPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=0 && z<lowPad[2] )
		return 25;
	// bottom front right: corresponds to top back left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=0 && y<lowPad[1] &&
		z>=0 && z<lowPad[2] )
		return 26;
	// bottom back right: corresponds to top front left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=0 && z<lowPad[2] )
		return 23;
	// bottom back left: corresponds to top front right
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=0 && z<lowPad[2] )
		return 24;
	// top front left: corresponds to bottom back right
	if( x>=0 && x<lowPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 21;
	// top front right: corresponds to bottom back left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=0 && y<lowPad[1] &&
		z>enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 22;
	// top back right: corresponds to bottom front left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 19;
	// top back left: corresponds to bottom front right
	if( x>=0 && x<lowPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 20;

	return 0;

}// end function


int
getRegionIndex( int regionType, int nReturnedRegion, float* receivedRegionID )
{
	int index = -1;

	for( int i=0; i<nReturnedRegion; i++ )
	{
		if( receivedRegionID[i] == (float)regionType )
		{
			index = i;
			break;
		}
	}

	return index;
}// end function

//
// makes MPI datatype for receiving one item
//
// cts: pointer to counts message
//
// side effects: allocates MPI datatype
//
// returns: pointer to MPI datatype
//
DIY_Datatype* RecvItemType(int *cts) {

  //MPI_Datatype *dtype = new MPI_Datatype;
  //MPI_Type_contiguous(4, MPI_FLOAT, dtype);

  //return dtype;
	DIY_Datatype *dtype = (DIY_Datatype *)malloc(sizeof(DIY_Datatype));

  struct map_block_t map[1] = {
    {DIY_FLOAT, OFST, 6*64*64+1, 0 },
  };
  DIY_Create_struct_datatype( 0, 1, map, dtype );

  return dtype;


}
//
// makes an MPI datatype for sending one item
//
// cts: pointer to counts message
// pts: pointer to points message
//
// side effects: allocates MPI datatype
//
// returns: pointer to MPI datatype
//
//
DIY_Datatype* SendItemType(int *cts, char** pts) {

  //MPI_Datatype *dtype = new MPI_Datatype; // datatype for one point
  //MPI_Type_contiguous(4, MPI_FLOAT, dtype);
  DIY_Datatype *dtype = (DIY_Datatype *)malloc(sizeof(DIY_Datatype));

  struct map_block_t map[1] = {
    {DIY_FLOAT, OFST, 6*64*64+1, 0},
  };
  DIY_Create_struct_datatype(0, 1, map, dtype);

  return dtype;
}

