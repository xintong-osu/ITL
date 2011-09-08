/**
 * @file MainIT_diy_local.cpp
 * Application program for local entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon 
 */

#include <mpi.h>

#include "diy.h"
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
#include "ITL_histogram.h"
#include "ITL_histogramconstants.h"
#include "ITL_field_regular.h"
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

ITL_field_regular<SCALAR> *enhancedScalarField = NULL;
//ITL_field_regular<VECTOR3> *vectorField = NULL;
ITL_field_regular<SCALAR> *scalarFieldBlockArray = NULL;

ITL_localentropy<SCALAR> *localEntropyComputer_scalar = NULL;
//ITL_localentropy<VECTOR3> *localEntropyComputer_vector = NULL;

int tot_blocks = 512;

double execTime[3];
clock_t starttime, endtime;

MPI_Datatype* SendItemType(int *cts, char** pts);
MPI_Datatype* RecvItemType(int *cts);
void Compute( int lid, int neighborhoodSize );
void Compute_oneblock( ITL_field_regular<SCALAR> *scalarFieldBlock,
					   float **dataToSend, int *nDataPointsToSend,
					   float **remotePoint, int neighborhoodSize );
void expandBlockData( float *blockData, float *lowF, float *highF,
					  float **receivedItems, float *enhancedData,
					  int neighborhoodSize );
int getRegionType( int x, int y, int z,
				   int *enhancedSize,
				   int neighborhoodSize );

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
	int given[3] = {0, 0, 0};
	int neighborhoodSize = 6;

	int fieldType = -1;
	int method = 0;
	int verboseMode = 1;

	float lowF[3];
	float highF[3];
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
	ITL_histogram::ITL_init_histogram( patchFile, nBin );

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
		if( fieldType == 0 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(data[i]));
		//if( fieldType == 1 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(vectordata[i]));

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
	DIY_Read_all();
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Allocate memory for storing received ghost layers and the local entropy fields for each block
	float *localEntropyList[nblocks];
	float *enhancedData[nblocks];
  	void **receivedItems; // received items from neighbors (generic pointer to a byte, not a string)

	// Compute ghost layers for all blocks
  	ITL_field_regular<SCALAR> scalarFieldBlockArray[nblocks];
	for( int k=0; k<nblocks; k++ )
	{
		// Compute block extent (without ghost layers)
		lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
		lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
		lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

		// Create a temporary block from original data (without ghost layers)
		scalarFieldBlockArray[k].setBounds( lowF, highF );
		scalarFieldBlockArray[k].setDataFull( data[k] );

		// Create chunks of data to be sent for this block
		Compute( k, neighborhoodSize );
	}

	// Scan through the blocks, compute local entropy fields
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<SCALAR> *scalarField = NULL;
	for( int k=0; k<nblocks; k++ )
	{
		if( fieldType == 0 )
		{

			// exchange neighbors
	   	 	DIY_Exchange_neighbors( &receivedItems, 1.0, &RecvItemType, &SendItemType );

			// Expand block data with ghost layers
			expandBlockData( data[k], lowF, highF, (float**)receivedItems,
							 enhancedData[k], neighborhoodSize );

   		 	// flush any remaining messages
  		  	// if multiple rounds of compute / exchange neighbors, call FlushNeighbors
  		  	// once for each time block, after the rounds complete
 		   	DIY_Flush_neighbors(&receivedItems, &RecvItemType );

			// Initialize ITL scalar field with block data
			scalarField = new ITL_field_regular<SCALAR>( enhancedData[k], nDim, lowF, highF );

			// Initialize class that can compute entropy
			localEntropyComputer_scalar = new ITL_localentropy<SCALAR>( scalarField );

			// Create bin field
			localEntropyComputer_scalar->computeHistogramBinField( "scalar", nBin );

			// Compute entropy
			localEntropyComputer_scalar->computeEntropyOfField( nBin, false );

			// Save local entropy field
			localEntropyList[k] = new float[scalarField->grid->nVertices];
			memcpy( localEntropyList[k], localEntropyComputer_scalar->getEntropyField(), sizeof(SCALAR)*scalarField->grid->nVertices );
			if( verboseMode == 1 ) printf( "Block Limits: %d, %d, %d, %d, %d, %d, local entropy computed.\n", scalarField->grid->lowInt[0],
															scalarField->grid->highInt[0],
															scalarField->grid->lowInt[1],
															scalarField->grid->highInt[1],
															scalarField->grid->lowInt[2],
															scalarField->grid->highInt[2] );

			// Clear up
			delete localEntropyComputer_scalar;
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
	execTime[1] = ITL_util<float>::endTimer( starttime );


	/*
	// Write local entropy field
	if( verboseMode == 1 ) printf( "Writing local entropy field ...\n" );
	starttime = ITL_util<float>::startTimer();
	MPI_Datatype *dtype = new MPI_Datatype; // datatype for output
	MPI_Type_contiguous( 1, MPI_FLOAT, dtype );
	//MPI_Type_contiguous( nblocks, MPI_FLOAT, dtype );
	MPI_Type_commit( dtype );
	MPI_Barrier( MPI_COMM_WORLD ); // everyone synchronizes again

	float **listptr = new float*[nblocks];
	for( int i=0; i<nblocks; i++ )
	{
		listptr[i] = localEntropyList[i];
	}

	DIY_Write_open_all( outFile, 0 );
	DIY_Write_blocks_all( (void **)&localEntropyList[0], nblocks, dtype );
	DIY_Write_close_all();
	execTime[2] = ITL_util<float>::endTimer( starttime );

	if( verboseMode == 1 ) printf( "%d: Read/Compute/Write Time: %f %f %f seconds\n", rank, execTime[0], execTime[1], execTime[2] );
	else printf( "%d, %f, %f, %f\n", rank, execTime[0], execTime[1], execTime[2] );
	*/

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
//
void Compute( int lid, int neighborhoodSize )
{
	int nItemToSendPerBlock = 26;
	float *dataToSend[nItemToSendPerBlock]; 
	int nDataPointsToSend[nItemToSendPerBlock]; 
	float *remotePoint[nItemToSendPerBlock];
	for( int i=0; i<nItemToSendPerBlock; i++ )
		remotePoint[nItemToSendPerBlock] = new float[3];

	// local computation here, producing an item that needs to be sent to 
  	// a neighboring block
   	Compute_oneblock( scalarFieldBlockArray+lid, dataToSend,
   					  nDataPointsToSend, &remotePoint[0], neighborhoodSize );
	
	// enqueue the item for sending to neighbor
	// in this example the point by which the neighbor is identified (5th arg.)
	// is the same as the item (2nd arg.), but this need not be the case
	// because the item can be any generic data
	  	
	for( int j=0; j<nItemToSendPerBlock; j++ )
	{
		DIY_Enqueue_item( lid, dataToSend[j], NULL, nDataPointsToSend[j] * sizeof(float), remotePoint[j] );
	}// end inner for

}// end for

// Here I assume x spans from left to right, y spans from bottom to top and z spans from front to back	
// First we cover 6 ghost faces (id:1 to 6)
// Next we cover 12 ghost edges (id:7 to 18)
// Then we cover 8 ghost corners (id:19 to 26)
void Compute_oneblock( ITL_field_regular<SCALAR> *scalarFieldBlock,
					   float **dataToSend, int *nDataPointsToSend,
					   float **remotePoint, int neighborhoodSize )
{
	int lowBoundary[3];
	int highBoundary[3];

	// 6 ghost faces
	// First we iterate through the side faces in a counterclockwise manner
	// Left (id: 1)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0] + neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[0] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[0], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[0] );
	remotePoint[0][0] = lowBoundary[0] - 1;
	remotePoint[0][1] = lowBoundary[1] + 1;
	remotePoint[0][2] = lowBoundary[2] + 1;
	
	// Front (id: 2)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[1] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[1], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[1] );
	remotePoint[1][0] = lowBoundary[0] + 1;
	remotePoint[1][1] = lowBoundary[1] + 1;
	remotePoint[1][2] = lowBoundary[2] - 1;

	// Right (id: 3)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0]-neighborhoodSize;;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[2] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[2], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[2] );
	remotePoint[2][0] = highBoundary[0] + 1;
	remotePoint[2][1] = highBoundary[1] - 1;
	remotePoint[2][2] = highBoundary[2] - 1;

	// Back (id: 4)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[3] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[3], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[3] );
	remotePoint[3][0] = highBoundary[0] - 1;
	remotePoint[3][1] = highBoundary[1] - 1;
	remotePoint[3][2] = highBoundary[2] + 1;

	// Now we do the bottom and the top face
	// Bottom (id: 5)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[4] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[4], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[4] );
	remotePoint[4][0] = highBoundary[0] - 1;
	remotePoint[4][1] = highBoundary[1] + 1;
	remotePoint[4][2] = highBoundary[2] - 1;


	// Top (id: 6)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[5] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[5], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[5] );
	remotePoint[5][0] = lowBoundary[0] + 1;
	remotePoint[5][1] = lowBoundary[1] - 1;
	remotePoint[5][2] = lowBoundary[2] + 1;
	
	// 12 ghost bars grazing the edges 
	// First we do the bottom 4 edges in a counterclockwise manner	
	// Bottom left (id: 7)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[6] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[6], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[6] );
	remotePoint[6][0] = lowBoundary[0] - 1;
	remotePoint[6][1] = lowBoundary[1] - 1;
	remotePoint[6][2] = lowBoundary[2] + 1;

	// Bottom front (id: 8)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[7] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[7], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[7] );
	remotePoint[7][0] = lowBoundary[0] + 1;
	remotePoint[7][1] = lowBoundary[1] - 1;
	remotePoint[7][2] = lowBoundary[2] - 1;

	// Bottom right (id: 9)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0]-neighborhoodSize;;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[8] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[8], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[8] );
	remotePoint[8][0] = highBoundary[0] + 1;
	remotePoint[8][1] = lowBoundary[1] - 1;
	remotePoint[8][2] = lowBoundary[2] + 1;

	// Bottom back (id: 10)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[9] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[9], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[9] );
	remotePoint[9][0] = lowBoundary[0] + 1;
	remotePoint[9][1] = lowBoundary[1] - 1;
	remotePoint[9][2] = highBoundary[2] + 1;

	// Next we do the top 4 edges in a counterclockwise manner	
	// Top left (id: 11)	
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
  	highBoundary[2] = scalarFieldBlock->grid->highInt[1];
	nDataPointsToSend[10] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[10], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[10] );
	remotePoint[10][0] = lowBoundary[0] - 1;
	remotePoint[10][1] = highBoundary[1] + 1;
	remotePoint[10][2] = lowBoundary[2] + 1;

	// Top front (id: 12)	
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[11] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[11], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[11] );
	remotePoint[11][0] = lowBoundary[0] + 1;
	remotePoint[11][1] = highBoundary[1] + 1;
	remotePoint[11][2] = lowBoundary[2] - 1;

	// Top right (id: 13)	
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[12] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[12], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[12] );
	remotePoint[12][0] = highBoundary[0] + 1;
	remotePoint[12][1] = highBoundary[1] + 1;
	remotePoint[12][2] = lowBoundary[2] + 1;

	// Top back (id: 14)	
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[13] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[13], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[13] );
	remotePoint[13][0] = lowBoundary[0] + 1;
	remotePoint[13][1] = highBoundary[1] + 1;
	remotePoint[13][2] = highBoundary[2] + 1;

	// Next we do the standing 4 edges in a counterclockwise manner	
	// Standing front left (id: 15) 
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
  	highBoundary[2] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	nDataPointsToSend[14] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[14], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[14] );
	remotePoint[14][0] = lowBoundary[0] - 1;
	remotePoint[14][1] = lowBoundary[1] + 1;
	remotePoint[14][2] = lowBoundary[2] - 1;

	// Standing front right (id: 16)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[15] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[15], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[15] );
	remotePoint[15][0] = highBoundary[0] + 1;
	remotePoint[15][1] = lowBoundary[1] + 1;
	remotePoint[15][2] = lowBoundary[2] - 1;

	// Standing back right (id: 17)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[16] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[16], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[16] );
	remotePoint[16][0] = highBoundary[0] + 1;
	remotePoint[16][1] = lowBoundary[1] + 1;
	remotePoint[16][2] = highBoundary[2] + 1;

	// Standing back left (id: 18)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[17] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[17], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[17] );
	remotePoint[17][0] = lowBoundary[0] - 1;
	remotePoint[17][1] = lowBoundary[1] + 1;
	remotePoint[17][2] = highBoundary[2] + 1;

	// 8 corners
	// First we do the bottom 4 corners in a counterclockwise manner	
	// Bottom front left (id: 19)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2]+neighborhoodSize;
	dataToSend[18] = scalarFieldBlock->getDataBetween( lowBoundary, highBoundary  );
	nDataPointsToSend[18] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[18], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[18] );
	remotePoint[18][0] = lowBoundary[0] - 1;
	remotePoint[18][1] = lowBoundary[1] - 1;
	remotePoint[18][2] = lowBoundary[2] - 1;

	// Bottom front right (id: 20)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[19] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[19], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[19] );
	remotePoint[19][0] = highBoundary[0] - 1;
	remotePoint[19][1] = highBoundary[1] + 1;
	remotePoint[19][2] = lowBoundary[2] - 1;

	// Bottom back right (id: 21)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[20] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[20], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[20] );
	remotePoint[20][0] = highBoundary[0] + 1;
	remotePoint[20][1] = lowBoundary[1] - 1;
	remotePoint[20][2] = highBoundary[2] + 1;

	// Bottom back left (id: 22)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->lowInt[1];
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;;
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->lowInt[1]+neighborhoodSize;;
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[21] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[21], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[21] );
	remotePoint[21][0] = lowBoundary[0] - 1;
	remotePoint[21][1] = lowBoundary[1] - 1;
	remotePoint[21][2] = highBoundary[2] + 1;

	// Next we do the top 4 corners in a counterclockwise manner	
	// Top front left (id: 23)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->lowInt[2]+neighborhoodSize;
	nDataPointsToSend[22] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[22], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[22] );
	remotePoint[22][0] = lowBoundary[0] - 1;
	remotePoint[22][1] = highBoundary[1] + 1;
	remotePoint[22][2] = lowBoundary[2] - 1;

	// Top front right (id: 24)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->lowInt[2];
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2]+neighborhoodSize;
	nDataPointsToSend[23] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[23], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[23] );
	remotePoint[23][0] = highBoundary[0] + 1;
	remotePoint[23][1] = highBoundary[1] + 1;
	remotePoint[23][2] = lowBoundary[2] - 1;


	// Top back right (id: 25)
	lowBoundary[0] = scalarFieldBlock->grid->highInt[0]-neighborhoodSize;
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->highInt[0];
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[24] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[24], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[24] );
	remotePoint[24][0] = highBoundary[0] + 1;
	remotePoint[24][1] = highBoundary[1] + 1;
	remotePoint[24][2] = highBoundary[2] + 1;

	// Top back left (id: 26)
	lowBoundary[0] = scalarFieldBlock->grid->lowInt[0];
  	lowBoundary[1] = scalarFieldBlock->grid->highInt[1]-neighborhoodSize;
	lowBoundary[2] = scalarFieldBlock->grid->highInt[2]-neighborhoodSize;
	highBoundary[0] = scalarFieldBlock->grid->lowInt[0]+neighborhoodSize;
  	highBoundary[1] = scalarFieldBlock->grid->highInt[1];
	highBoundary[2] = scalarFieldBlock->grid->highInt[2];
	nDataPointsToSend[25] = ( highBoundary[0] - lowBoundary[0] ) *
						   ( highBoundary[1] - lowBoundary[1] ) *
						   ( highBoundary[2] - lowBoundary[2] );
	memcpy( dataToSend[25], scalarFieldBlock->getDataBetween( lowBoundary, highBoundary ),
			sizeof(SCALAR)*nDataPointsToSend[25] );
	remotePoint[25][0] = lowBoundary[0] - 1;
	remotePoint[25][1] = highBoundary[1] + 1;
	remotePoint[25][2] = highBoundary[2] + 1;

}// end for

// Here I assume x spans from left to right, y spans from bottom to top and z spans from front to back
void expandBlockData( float *blockData, float *lowF, float *highF,
					  float **receivedItems, float *enhancedData,
					  int neighborhoodSize )
{
	int enhancedSize[3];
	int nItemRecvdPerBlock = 26;
	int recvdItemIndex[nItemRecvdPerBlock];
	for( int i=0; i<nItemRecvdPerBlock; i++ )
		recvdItemIndex[i] = 0;
	
	int indexid = 0;
	int ownDataIndex = 0;
	for( int z=0; z<enhancedSize[2]; z++ )
	{	
		for( int y=0; y<enhancedSize[1]; y++ )
		{
			for( int x=0; x<enhancedSize[0]; x++ )
			{
				int regionType = getRegionType( x, y, z, enhancedSize, neighborhoodSize );

				if( regionType == 0 )
				{
					enhancedData[indexid] = blockData[ownDataIndex];
					ownDataIndex++;
				}
				else
				{
					enhancedData[indexid] = receivedItems[regionType-1][recvdItemIndex[regionType-1]];
					recvdItemIndex[regionType-1]++;
				}


				indexid++;
 
			}// end for z
		}// end for y
	}// end for z				
	
}// end function 

// Here I assume x spans from left to right, y spans from bottom to top and z spans from front to back
int getRegionType( int x, int y, int z, int *enhancedSize, int neighborhoodSize ) 
{
	// Determine if point in any of the 6 ghost faces
	// left	: corresponds to right
	if( x<=0 && x<neighborhoodSize &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 3;
	// front: corresponds to back
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=0 && z<enhancedSize[2]-neighborhoodSize ) return 4;
	// right: corresponds to left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 1;
	// back: corresponds to front
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 2;
	// bottom: corresponds to top
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=0 && y<enhancedSize[1]-neighborhoodSize &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 6;
	// top: corresponds to bottom
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 5;
	 
	// Determine if point in any of the 12 ghost edge-grazing bars
	// bottom left	: corresponds to top right
	if( x>=0 && x<neighborhoodSize &&
		y>=0 && y<neighborhoodSize &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 13;
	// bottom front: corresponds to top back
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=0 && y<neighborhoodSize &&
		z>=0 && z<neighborhoodSize ) return 14;
	// bottom right: corresponds to top left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=0 && y<enhancedSize[1]-neighborhoodSize &&
		z>=neighborhoodSize && z<neighborhoodSize ) return 11;
	// bottom back: corresponds to top front
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=0 && y<enhancedSize[1]-neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 12;
	// top left	: corresponds to bottom right
	if( x>=0 && x<neighborhoodSize &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=neighborhoodSize && z<enhancedSize[2]-neighborhoodSize ) return 9;
	// top front: corresponds to bottom back
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=0 && z<neighborhoodSize ) return 10;
	// top right: corresponds to bottom left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=neighborhoodSize && z<neighborhoodSize ) return 7;
	// top back: corresponds to bottom front
	if( x>=neighborhoodSize && x<enhancedSize[0]-neighborhoodSize &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 8;
	// standing front left: corresponds to standing back right
	if( x>=0 && x<neighborhoodSize &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=0 && z<neighborhoodSize ) return 17;
	// standing front right: corresponds to standing back left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=0 && z<neighborhoodSize ) return 18;
	// standing back right: corresponds to front left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 15;
	// standing back left: corresponds to front front right
	if( x>=0 && x<neighborhoodSize &&
		y>=neighborhoodSize && y<enhancedSize[1]-neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 16;

	// Determine if point in any of the 8 ghost corners		
	// bottom front left: corresponds to top back right
	if( x>=0 && x<neighborhoodSize &&
		y>=0 && y<neighborhoodSize &&
		z>=0 && z<neighborhoodSize ) return 25;
	// bottom front right: corresponds to top back left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=0 && y<neighborhoodSize &&
		z>=0 && z<neighborhoodSize ) return 26;
	// bottom back right: corresponds to top front left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=0 && y<neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 23;
	// bottom back left: corresponds to top front right
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=0 && y<neighborhoodSize &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 24;
	// top front left: corresponds to bottom back right
	if( x>=0 && x<neighborhoodSize &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=0 && z<neighborhoodSize ) return 21;
	// top front right: corresponds to bottom back left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=0 && z<neighborhoodSize ) return 22;
	// top back right: corresponds to bottom front left
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 19;
	// top back left: corresponds to bottom front right
	if( x>=enhancedSize[0]-neighborhoodSize && x<enhancedSize[0] &&
		y>=enhancedSize[1]-neighborhoodSize && y<enhancedSize[1] &&
		z>=enhancedSize[2]-neighborhoodSize && z<enhancedSize[2] ) return 20;

	return 0;

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
MPI_Datatype* RecvItemType(int *cts) {

  MPI_Datatype *dtype = new MPI_Datatype;
  MPI_Type_contiguous(4, MPI_FLOAT, dtype);

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
MPI_Datatype* SendItemType(int *cts, char** pts) {

  MPI_Datatype *dtype = new MPI_Datatype; // datatype for one point
  MPI_Type_contiguous(4, MPI_FLOAT, dtype);

  return dtype;

}

