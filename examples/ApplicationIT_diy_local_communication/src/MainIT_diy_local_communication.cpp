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

int dataSize[3];
int blockSize[3];
int tot_blocks = 512;
int nblocks;				// My local number of blocks

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR> *histMapper_scalar = NULL;

ITL_field_regular<int> *binField = NULL;
ITL_field_regular<SCALAR> *enhancedScalarField = NULL;
ITL_field_regular<SCALAR>** scalarFieldBlockArray = NULL;
ITL_field_regular<VECTOR3>** vectorFieldBlockArray = NULL;

ITL_localentropy<SCALAR>* localEntropyComputer_scalar = NULL;
ITL_localentropy<VECTOR3>* localEntropyComputer_vector = NULL;

// Nearest neighbor exchange related variables
int neighborhoodSize = 6;
int maxDimLength = 64;
int nItemToSendPerBlock = 26;
int maxGhostSize;
SCALAR*** dataToSend = NULL;
//SCALAR* dataToSend = NULL;
//void ***receivedItems = NULL;
unsigned char dirs[26] = { DIY_X0, DIY_Y0, DIY_X1, DIY_Y1, DIY_Z0, DIY_Z1,
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

double execTime[4];
clock_t starttime, endtime;

int verboseMode = 1;

DIY_Datatype* SendItemType(int *cts, char** pts);
DIY_Datatype* RecvItemType(int *cts);
void* CreateWriteType( void *item, int lid, DIY_Datatype *dtype );
void Compute( int did, int lid );
void expandBlockData( float **blockData,
					  float *lowF, float *highF,
		  	  	  	  int* lowPad, int* highPad,
		  	  	  	  int* paddedLow, int* paddedHigh,
		  	  	  	  float **receivedItems,
		  	  	  	  int nReceivedItems, //float* receivedItemId,
		  	  	  	  float **enhancedData,
		  	  	  	  int neighborhoodSize );
int getRegionType( int x, int y, int z,
				   int *enhancedSize,
				   int* lowPad, int* highPad);
int getRegionIndex( int regionType,
					int nReturnedRegion, float* receivedRegionID );

/**
 * Main function.
 * Program starts from here.
 */
int main( int argc, char** argv )
{
	int numProcs;
	int rank;
	int num_threads = 4; // number of threads DIY can use

	int nDim;

	int nBin;
	float histogramLowEnd = 0;
	float histogramHighEnd = 0;

	float lowF[3];
	float highF[3];
	int given[3] = {0, 0, 0};
	int lowPad[3] = {0, 0, 0};
	int highPad[3] = {0, 0, 0};
	int paddedLow[3] = {0, 0, 0};
	int paddedHigh[3] = {0, 0, 0};
	int neighborhoodSizeArray[3] = {0, 0, 0};

	int fieldType = -1;
	int method = 0;

	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

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
	neighborhoodSize = atoi( ITL_util<float>::getArgWithName( "neighborhoodSize", &argNames, &argValues ) );
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
	DIY_Init( nDim, dataSize, num_threads, MPI_COMM_WORLD );
	if( verboseMode == 1 )	fprintf( stderr, "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Decompose domain (with ghost layers)
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int did = DIY_Decompose( ROUND_ROBIN_ORDER, tot_blocks, &nblocks, 0, 0, given );

	// Allocate memory for pointers that will hold block data
	SCALAR* data[nblocks];
	VECTOR3* vectordata[nblocks];
	if( fieldType == 0 ) memset( data, 0, sizeof(float*) * nblocks );
	//if( fieldType == 1 ) memset( vectordata, 0, sizeof(VECTOR3*) * nblocks );

	// Serially visit the blocks (?)
	int* diy_min = new int[3*nblocks];
	int* diy_max = new int[3*nblocks];
	int* diy_size = new int[3*nblocks];

	starttime = ITL_util<float>::startTimer();
	int maxDimLength = 0;
	for (int i = 0; i < nblocks; i++)
	{
		DIY_Block_starts_sizes( did, i, &diy_min[3*i], &diy_size[3*i] );

		// post a read for the block
		if( fieldType == 0 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, DIY_FLOAT, (void**)&(data[i]));
		//if( fieldType == 1 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(vectordata[i]));

		// print the block bounds
		for (int j = 0; j < 3; j++)
			diy_max[3*i+j] = diy_min[3*i+j] + diy_size[3*i+j] - 1;

		if( verboseMode == 1 )
			fprintf( stderr, "process rank = %d "
					 "block local id = %d "
					 "min = [%d %d %d] "
					 "max = [%d %d %d] "
					 "size = [%d %d %d]\n",
					 rank, i,
					 diy_min[3*i], diy_min[3*i+1], diy_min[3*i+2],
					 diy_max[3*i], diy_max[3*i+1], diy_max[3*i+2],
					 diy_size[3*i], diy_size[3*i+1], diy_size[3*i+2] );

		// Compute largest face size for this block
		if( diy_size[3*i] > maxDimLength ) maxDimLength =  diy_size[3*i];
		if( diy_size[3*i+1] > maxDimLength ) maxDimLength =  diy_size[3*i+1];
		if( diy_size[3*i+2] > maxDimLength ) maxDimLength =  diy_size[3*i+2];
	}
	if( verboseMode == 1 )	fprintf( stderr, "Locally computed largest face size: %d\n", maxDimLength );

	// Read actual data (everyone synchronizes after reading data)
	DIY_Read_data_all();
	if( verboseMode == 1 )	fprintf( stderr, "Data read complete .. all processes synced\n" );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Determine global largest face size
	int* sendMaxLocalBlockSize = new int[numProcs];
	int* recvMaxLocalBlockSize = new int[numProcs];
	memset( sendMaxLocalBlockSize, 0, numProcs*sizeof(int) );
	sendMaxLocalBlockSize[rank] = maxDimLength;
	// Start time
	starttime = ITL_util<float>::startTimer();
	MPI_Allreduce( sendMaxLocalBlockSize, recvMaxLocalBlockSize, numProcs, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    maxDimLength = recvMaxLocalBlockSize[rank];
    //if( verboseMode == 1 )
    if( verboseMode == 1 )	fprintf( stderr, "Globally computed largest face size: %d\n", maxDimLength );
    MPI_Barrier( MPI_COMM_WORLD );
    execTime[2] = ITL_util<float>::endTimer( starttime );
    // Stop time
    maxGhostSize = neighborhoodSize * maxDimLength * maxDimLength;
    delete [] sendMaxLocalBlockSize;
    delete [] recvMaxLocalBlockSize;

	// Compute ghost layers for all blocks
	scalarFieldBlockArray = new ITL_field_regular<SCALAR>*[nblocks];
	if( verboseMode == 1 )	fprintf( stderr, "Creating local fields ...\n" );
	for( int k=0; k<nblocks; k++ )
	{
		// Compute block extent (without ghost layers)
		lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
		lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
		lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

		#ifdef DEBUG_MODE
		float m = ITL_util<float>::Min( data[k], (highF[0] - lowF[0] + 1)*(highF[1] - lowF[1] + 1)*(highF[2] - lowF[2] + 1) );
		float M = ITL_util<float>::Max( data[k], (highF[0] - lowF[0] + 1)*(highF[1] - lowF[1] + 1)*(highF[2] - lowF[2] + 1) );
		fprintf( stderr, "Block value range (no ghost): %g %g\n", m, M );
		#endif

		// Create a temporary block from original data (without ghost layers)
		scalarFieldBlockArray[k] = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );
	}
	if( verboseMode == 1 ) fprintf( stderr, "Creating local fields done ...\n" );

	// Allocate memory for packets to send
	//dataToSend = new SCALAR[ nblocks * nItemToSendPerBlock * (4 + maxGhostSize) ];

	starttime = ITL_util<float>::startTimer();
	dataToSend = new SCALAR**[nblocks];
	for( int k=0; k<nblocks; k++ )
	{
		dataToSend[k] = new SCALAR*[nItemToSendPerBlock];
		for( int i=0; i<nItemToSendPerBlock; i++ )
		   dataToSend[k][i] = new SCALAR[4 + maxGhostSize];
	}

	// Create chunks of data to be sent
	for( int k=0; k<nblocks; k++ )
	{
		if( verboseMode == 1 ) fprintf( stderr, "Entering compute for block %d...\n", k );
		Compute( rank, k );
	}
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Scan through the blocks, compute local entropy fields
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<SCALAR> *scalarField = NULL;

	// Allocate memory for storing received ghost layers and

	SCALAR* enhancedData[nblocks];
	memset( enhancedData, 0, sizeof(SCALAR*) * nblocks );
	void*** receivedItems = new void**[nblocks]; // received items from neighbors (generic pointer to a byte, not a string)
	int* num_items_recvd = new int[nblocks]; // number of received items in each block
	float *returnedItemId = new float[nblocks]; // number of received items in each block

	// exchange neighbors
	if( verboseMode == 1 ) fprintf( stderr, "Starting Neighbor exchange ...\n" );
  	starttime = ITL_util<float>::startTimer();
  	//DIY_Exchange_neighbors( did, receivedItems, num_items_recvd, 1.0, &RecvItemType, &SendItemType );
  	execTime[2] = execTime[2] + ITL_util<float>::endTimer( starttime );
  	if( verboseMode == 1 ) fprintf( stderr, "Neighbor exchange done ...\n" );

	// the local entropy fields for each block
	// The memset to 0 is needed to tell DIY to allocate the memory for us
	SCALAR* localEntropyList[nblocks];
	memset( localEntropyList, 0, sizeof(SCALAR*) * nblocks );
	int enhancedSize[3];
	int nVertEnhanced;
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

				enhancedSize[i] = ( paddedHigh[i] - paddedLow[i] + 1 );

			}// end for

			#ifdef DEBUG_MODE
			fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																	highF[0], highF[1], highF[2] );
			fprintf( stderr, "Padded Block boundary: %d %d %d %d %d %d\n", paddedLow[0], paddedLow[1], paddedLow[2],
																		   paddedHigh[0], paddedHigh[1], paddedHigh[2] );
			fprintf( stderr, "Number of items received by block %d: %d \n", k, num_items_recvd[k] );
			#endif

			// Expand block data with ghost layers
			if( verboseMode == 1 ) fprintf( stderr, "expanding block data for block %d ...\n", k );
			expandBlockData( &data[k],
							 lowF, highF,
							 lowPad, highPad,
							 paddedLow,paddedHigh,
							 (float**)receivedItems[k],
							 num_items_recvd[k],// returnedItemId,
							 &enhancedData[k], neighborhoodSize );

			// Initialize ITL scalar field with block data
			if( verboseMode == 1 ) fprintf( stderr, "creating scalar field for block %d ...\n", k );
			scalarField = new ITL_field_regular<SCALAR>( enhancedData[k], nDim, lowF, highF,
														 lowPad, highPad, neighborhoodSizeArray );

			#ifdef DEBUG_MODE
			fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																	highF[0], highF[1], highF[2] );
			fprintf( stderr, "Pad: %d %d %d %d %d %d\n", lowPad[0], lowPad[1], lowPad[2],
														highPad[0], highPad[1], highPad[2] );
			fprintf( stderr, "Padded size: %d %d %d\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );
			fprintf( stderr, "Nhood: %d %d %d\n", neighborhoodSizeArray[0], neighborhoodSizeArray[1], neighborhoodSizeArray[2] );
			#endif

			// Create bin field
			if( verboseMode == 1 ) fprintf( stderr, "Creating histogram mapper for block %d ...\n", k );
			if( k == 0 && histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBin );

			#ifdef DEBUG_MODE
			nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];
			FILE* debugFile1;
			stringstream temp1( stringstream::in | stringstream::out );
			temp1 << "./debug/block_" << rank << "_" << k << ".bin";
			string filename1 = temp1.str();
			debugFile1 = fopen( filename1.c_str(), "wb" );
			assert( debugFile1 != NULL );
			fwrite( enhancedSize, sizeof(int), 3, debugFile1 );
			fwrite( enhancedData[k], sizeof(float), nVertEnhanced, debugFile1 );
			fclose( debugFile1 );
			FILE* debugFile2;
			stringstream temp2( stringstream::in | stringstream::out );
			temp2 << "./debug/bblock_" << rank << "_" << k << ".bin";
			string filename2 = temp2.str();
			debugFile2 = fopen( filename2.c_str(), "wb" );
			nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];
			assert( debugFile2 != NULL );
			//fwrite( enhancedSize, sizeof(int), 3, dataFile );
			fwrite( binField->getDataFull(), sizeof(int), nVertEnhanced, debugFile2 );
			fclose( debugFile2 );

			int binm = ITL_util<int>::Min( binField->getDataFull(), nVertEnhanced  );
			int binM = ITL_util<int>::Max( binField->getDataFull(), nVertEnhanced );
			fprintf( stderr, "binField m M: %d %d\n", binm, binM );
			#endif

			// Initialize class that can compute entropy
			if( verboseMode == 1 ) fprintf( stderr, "creating entropy class for block %d ...\n", k );
			localEntropyComputer_scalar = new ITL_localentropy<SCALAR>( binField, histogram, nBin );

			// Compute entropy
			if( verboseMode == 1 ) fprintf( stderr, "computing entropy for block %d ...\n", k );
			localEntropyComputer_scalar->computeLocalEntropyOfField( false );

			// Save local entropy field
			if( verboseMode == 1 ) fprintf( stderr, "saving entropy field for block %d ...\n", k );
			#ifdef DEBUG_MODE
			fprintf( stderr, "outblock boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
																	highF[0], highF[1], highF[2] );
			fprintf( stderr, "Write buffer size for block %d: %d\n", k, blockSize[0]*blockSize[1]*blockSize[2]+6 );
			#endif
			blockSize[0] = (int)( highF[0] - lowF[0] + 1 );
			blockSize[1] = (int)( highF[1] - lowF[1] + 1 );
			blockSize[2] = (int)( highF[2] - lowF[2] + 1 );
			localEntropyList[k] = new SCALAR[blockSize[0]*blockSize[1]*blockSize[2]+6];

			localEntropyList[k][0] = diy_min[3*k];
			localEntropyList[k][1] = diy_min[3*k+1];
			localEntropyList[k][2] = diy_min[3*k+2];
			localEntropyList[k][3] = diy_max[3*k];
			localEntropyList[k][4] = diy_max[3*k+1];
			localEntropyList[k][5] = diy_max[3*k+2];
			memcpy( localEntropyList[k]+6,
					localEntropyComputer_scalar->getEntropyField()->getDataFull(),
					sizeof(SCALAR)*(blockSize[0]*blockSize[1]*blockSize[2]) );

			#ifdef DEBUG_MODE
			SCALAR em = ITL_util<float>::Min( localEntropyComputer_scalar->getEntropyField()->getDataFull(), blockSize[0]*blockSize[1]*blockSize[2]  );
			SCALAR eM = ITL_util<float>::Max( localEntropyComputer_scalar->getEntropyField()->getDataFull(), blockSize[0]*blockSize[1]*blockSize[2] );
			fprintf( stderr, "enField m M: %g %g\n", em, eM );
			#endif

			#ifdef DEBUG_MODE
			FILE* debugFile3;
			stringstream temp3( stringstream::in | stringstream::out );
			temp3 << "./debug/lblock_" << rank << "_" << k << ".bin";
			string filename3 = temp3.str();
			debugFile3 = fopen( filename3.c_str(), "wb" );
			assert( debugFile3 != NULL );
			fwrite( blockSize, sizeof(int), 3, debugFile3 );
			fwrite( &localEntropyList[k][6], sizeof(float), blockSize[0]*blockSize[1]*blockSize[2], debugFile3 );
			fclose( debugFile3 );
			#endif

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
	//fprintf( stderr, "Flushing neighbors ...\n" );
	//DIY_Flush_neighbors( receivedItems, num_items, &RecvItemType );
	//fprintf( stderr, "Done\n" );
	execTime[1] = execTime[1] + ITL_util<float>::endTimer( starttime );

	// Write local entropy field
	if( verboseMode == 1 ) printf( "Writing local entropy field ...\n" );
	starttime = ITL_util<float>::startTimer();
	//DIY_Write_open_all( outFile, 0 );
	//DIY_Write_blocks_all( (void **)&localEntropyList[0], nblocks, NULL, 0, &CreateWriteType );
	//DIY_Write_close_all();
	execTime[3] = ITL_util<float>::endTimer( starttime );


	if( verboseMode == 1 )
		fprintf( stderr, "%d: Read/Compute/Communicaiton/Write Time: %g %g %g %g seconds\n",
								   rank, execTime[0], execTime[1], execTime[2], execTime[3] );
	else
		fprintf( stderr, "%d, %g, %g, %g, %g\n", rank, execTime[0], execTime[1], execTime[2], execTime[3]  );

	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	for( int k = 0; k<nblocks; k++ )
		delete [] localEntropyList[k];
//	delete [] dataToSend;

	// Finalize MPI
	//if( fieldType == 1 ) MPI_Type_free( &complex );
	MPI_Finalize();

}// End main




//
// computation for my local block number lid (local ID)
// Here I assume x spans from left to right, y spans from front to back and z spans from bottom to top
//
void Compute( int did, int lid )
{
	float m, M;
	int lowBoundary[3];
	int highBoundary[3];
	int nDataPointsToSend[nItemToSendPerBlock];

	// local computation here, producing an item that needs to be sent to 
  	// a neighboring block
	if( verboseMode == 1 ) fprintf( stderr, "Entering compute one block for block %d ..\n", lid );

	for( int iN = 0; iN<nItemToSendPerBlock; iN++ )
	{
	    scalarFieldBlockArray[lid]->getBounds( lowBoundary, highBoundary );

		// 6 ghost faces
		// First we iterate through the side faces in a counterclockwise manner
		if( iN == 0 )
		{
			// Left (id: 1)
			highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		}
		else if( iN == 1 )
		{
			// Front (id: 2)
			highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
		}
		else if( iN == 2 )
		{
			// Right (id: 3)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		}
		else if( iN == 3 )
		{
			// Back (id: 4)
			lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		}
		else if( iN == 4 )
		{
			// Bottom (id: 5)
			highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 5 )
		{
			// Top (id: 6)
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		// 12 ghost bars grazing the edges
		// First we do the bottom 4 edges in a counterclockwise manner
		else if( iN == 6 )
		{
			// Bottom left (id: 7)
			highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 7 )
		{
			// Bottom front (id: 8)
		  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
			highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 8 )
		{
			// Bottom right (id: 9)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		   	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 9 )
		{
			// Bottom back (id: 10)
			lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		// Next we do the top 4 edges in a counterclockwise manner
		else if( iN == 10 )
		{
			// Top left (id: 11)
		    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		    lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 11 )
		{
			// Top front (id: 12)
			lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 12 )
		{
			// Top right (id: 13)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		  	lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 13 )
		{
			// Top back (id: 14)
		  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		// Next we do the standing 4 edges in a counterclockwise manner
		else if( iN == 14 )
		{
			// Standing front left (id: 15)
			highBoundary[0] = lowBoundary[0] + neighborhoodSize -1;
		  	highBoundary[1] = lowBoundary[1] + neighborhoodSize -1;
		}
		else if( iN == 15 )
		{
			// Standing front right (id: 16)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
			highBoundary[1] = lowBoundary[1] + neighborhoodSize -1;
		}
		else if( iN == 16 )
		{
			// Standing back right (id: 17)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		 	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		}
		else if( iN == 17 )
		{
			// Standing back left (id: 18)
		    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		    lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		}
		// 8 corners
		// First we do the bottom 4 corners in a counterclockwise manner
		else if( iN == 18 )
		{
			// Bottom front left (id: 19)
			highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
			highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 19 )
		{
			// Bottom front right (id: 20)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		  	highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
			highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 20 )
		{
			// Bottom back right (id: 21)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		 	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		else if( iN == 21 )
		{
			// Bottom back left (id: 22)
			highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
		  	highBoundary[2] = lowBoundary[2] + neighborhoodSize - 1;
		}
		// Next we do the top 4 corners in a counterclockwise manner
		else if( iN == 22 )
		{
			// Top front left (id: 23)
		    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		    highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
		    lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 23 )
		{
			// Top front right (id: 24)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
			highBoundary[1] = lowBoundary[1] + neighborhoodSize - 1;
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 24 )
		{
			// Top back right (id: 25)
			lowBoundary[0] = highBoundary[0] - neighborhoodSize + 1;
		  	lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}
		else if( iN == 25 )
		{
			// Top back left (id: 26)
		    highBoundary[0] = lowBoundary[0] + neighborhoodSize - 1;
		    lowBoundary[1] = highBoundary[1] - neighborhoodSize + 1;
			lowBoundary[2] = highBoundary[2] - neighborhoodSize + 1;
		}

		nDataPointsToSend[iN] = ( highBoundary[0] - lowBoundary[0] + 1 ) *
								( highBoundary[1] - lowBoundary[1] + 1 ) *
								( highBoundary[2] - lowBoundary[2] + 1 );
		#ifdef DEBUG_MODE
		fprintf( stderr, "%d-th boundary: %d %d %d %d %d %d\n", iN, lowBoundary[0], lowBoundary[1], lowBoundary[2],
																highBoundary[0], highBoundary[1], highBoundary[2] );
		fprintf( stderr, "%d-th pad size %d\n", iN, nDataPointsToSend[iN] );
		#endif

		dataToSend[lid][iN][0] = iN+1;						// Layer type
		dataToSend[lid][iN][1] = did;						// Sender rank
		dataToSend[lid][iN][2] = lid;						// Sender local block id
		dataToSend[lid][iN][3] = nDataPointsToSend[iN];		// Message size
		scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[lid][iN][4] );

		/*
		int startIndex = lid * nItemToSendPerBlock * (4+maxGhostSize) +
						 iN * (4+maxGhostSize);

		dataToSend[startIndex] = iN+1;							// Layer type
		dataToSend[startIndex+1] = rank;						// Sender rank
		dataToSend[startIndex+2] = lid;							// Sender local block id
		dataToSend[startIndex+3] = nDataPointsToSend[iN];		// Message size
		scalarFieldBlockArray[lid]->getDataBetween( lowBoundary, highBoundary, &dataToSend[startIndex+4] );
		*/

		#ifdef DEBUG_MODE
		m= ITL_util<float>::Min( &dataToSend[startIndex+4], nDataPointsToSend[iN] );
		M = ITL_util<float>::Max( &dataToSend[startIndex+4], nDataPointsToSend[iN] );
		fprintf( stderr, "%d-th packet size and range: %d %g %g\n", iN, nDataPointsToSend[iN], m, M );
		fprintf( stderr, "startindex: %d\n", startIndex );
		#endif

		DIY_Enqueue_item_dirs( did, lid, (void *)( &dataToSend[lid][iN][0] ), NULL,
							   (nDataPointsToSend[iN]+4) * sizeof(SCALAR),
							   &dirs[iN], 1, NULL );
		//DIY_Enqueue_item_dirs( lid, (void *)( &dataToSend[startIndex] ), NULL,
		//					   (nDataPointsToSend[iN]+4) * sizeof(SCALAR),
		//					   &dirs[iN], 1, NULL );

		if( verboseMode == 1 ) fprintf( stderr, "Item %d complete\n", iN );

	}// end if : traverse though ghost layers

	if( verboseMode == 1 ) fprintf( stderr, "Out from compute one block for block %d ..\n", lid );

}// end for

// Here I assume x spans from left to right, y spans from bottom to top and z spans from front to back
void expandBlockData( float **blockData,
					  float *lowF, float *highF,
					  int* lowPad, int* highPad,
					  int* paddedLow, int* paddedHigh,
					  float **receivedItems,
					  int nReceivedItems, //float* receivedItemId,
					  float **enhancedData,
					  int neighborhoodSize )
{
	float m, M;
	int blocksize[3], enhancedSize[3];
	int recvdItemIndexArray[nReceivedItems];
	float recvdItemId[nReceivedItems];

	for( int i=0; i<nReceivedItems; i++ )
		recvdItemIndexArray[i] = 4;
	
	for( int i=0; i<nReceivedItems; i++ )
	{
		recvdItemId[i] = ((float*)receivedItems[i])[0];
		#ifdef DEBUG_MODE
		fprintf( stderr, "%d-th recvd item is of type %g\n", i, recvdItemId[i] );
		#endif
	}

	// Allocate memory for enhanced block
	for( int i=0; i<3; i++ )
	{
		blocksize[i] = (int)( highF[i] - lowF[i] + 1 );
		enhancedSize[i] = ( paddedHigh[i] - paddedLow[i] + 1 );
	}
	int nVertEnhanced = enhancedSize[0] * enhancedSize[1] * enhancedSize[2];

	(*enhancedData) = new float[nVertEnhanced];

	#ifdef DEBUG_MODE
	fprintf( stderr, "Block boundary: %g %g %g %g %g %g\n", lowF[0], lowF[1], lowF[2],
															highF[0], highF[1], highF[2] );
	fprintf( stderr, "Pad size: %d %d %d %d %d %d\n", lowPad[0], lowPad[1], lowPad[2],
													  highPad[0], highPad[1], highPad[2] );
	fprintf( stderr, "Block size: %d %d %d\n", blockSize[0], blockSize[1], blockSize[2] );
	fprintf( stderr, "Enhanced size: %d %d %d\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );
	fprintf( stderr, "Actual block size: %ld\n", blockSize[0]*blockSize[1]*blockSize[2] );
	fprintf( stderr, "Enhanced block size: %ld\n", nVertEnhanced );
	#endif

	int indexid = 0, regionType = 0, itemIndex = 0;
	int ownDataIndex = 0;
	int mirrorDataCount = 0;
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
				}
				else
				// Inside one of the ghost layers
				{
					// Get Region Type
					regionType = getRegionType( x, y, z, enhancedSize, lowPad, highPad );

					// Get index of region
					itemIndex = getRegionIndex( regionType, nReceivedItems, recvdItemId );

					#ifdef DEBUG_MODE
					if( itemIndex < 0 )
					{
						fprintf( stderr, "block bounds: <%g %g %g><%g %g %g>\n", lowF[0], lowF[1], lowF[2], highF[0], highF[1], highF[2] );
						fprintf( stderr, "block bounds: <%d %d %d><%d %d %d>\n", lowPad[0], lowPad[1], lowPad[2],
																				 highPad[0], highPad[1], highPad[2] );
						fprintf( stderr, "enhanced size: <%d %d %d>\n", enhancedSize[0], enhancedSize[1], enhancedSize[2] );
						fprintf( stderr, "num recvd item: %d\n", nReceivedItems );
						fprintf( stderr, "%d %d %d\n", x, y, z );
						fprintf( stderr, "region type and item index: %d %d\n", itemIndex, regionType );
						for( int i=0; i<nReceivedItems; i++ )
							fprintf( stderr, "%g\n", recvdItemId[i] );
					}
					#endif

					if( itemIndex == -1 )
						mirrorDataCount ++;
					else
					{
						(*enhancedData)[indexid] = receivedItems[itemIndex][recvdItemIndexArray[itemIndex]];
						recvdItemIndexArray[itemIndex]++;
					}
				}

				indexid++;
 
			}// end for z
		}// end for y
	}// end for z
	
	#ifdef DEBUG_MODE
	int sum = 0;
	for( int i=0; i<nReceivedItems; i++ )
	{
		fprintf( stderr, "Ghost layer id %d, type: %g,  Points came: %g, from: %g/%g, Points read: %d\n",
				         i,
				         receivedItems[i][0], receivedItems[i][3],
				         receivedItems[i][1], receivedItems[i][2],
				         recvdItemIndexArray[i]-4 );
		sum += (recvdItemIndexArray[i]-4);
	}

	fprintf( stderr, "inside points: %d\n", ownDataIndex );
	fprintf( stderr, "mirror points: %d\n", mirrorDataCount );
	fprintf( stderr, "Total points written: %d\n", sum + ownDataIndex + mirrorDataCount );
	m = ITL_util<float>::Min( (*blockData), blocksize[0]*blocksize[1]*blocksize[2] );
	M = ITL_util<float>::Max( (*blockData), blocksize[0]*blocksize[1]*blocksize[2] );
	fprintf( stderr, "Block data range: %g %g\n", m, M );
	m = ITL_util<float>::Min( (*enhancedData), enhancedSize[0]*enhancedSize[1]*enhancedSize[2] );
	M = ITL_util<float>::Max( (*enhancedData), enhancedSize[0]*enhancedSize[1]*enhancedSize[2] );
	fprintf( stderr, "xm xM: %g %g\n", m, M );
	for( int i=0; i<nReceivedItems; i++ )
	{
		m = ITL_util<float>::Min( &receivedItems[i][4], recvdItemIndexArray[i]-4 );
		M = ITL_util<float>::Max( &receivedItems[i][4], recvdItemIndexArray[i]-4 );
		fprintf( stderr, "%d: rm rM: %g %g\n", i, m, M );
	}
	#endif


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
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=lowPad[2] && z<enhancedSize[2]-highPad[2] )
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
		return 13;
	// bottom front: corresponds to top back
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=0 && z<lowPad[2] )
		return 14;
	// bottom right: corresponds to top left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=0 && z<lowPad[2] )
		return 11;
	// bottom back: corresponds to top front
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=0 && z<lowPad[2] )
		return 12;
	/////////////////////////////////////////
	// top left	: corresponds to bottom right
	if( x>=0 && x<lowPad[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 9;
	// top front: corresponds to bottom back
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=0 && y<lowPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 10;
	// top right: corresponds to bottom left
	if( x>=enhancedSize[0]-highPad[0] && x<enhancedSize[0] &&
		y>=lowPad[1] && y<enhancedSize[1]-highPad[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 7;
	// top back: corresponds to bottom front
	if( x>=lowPad[0] && x<enhancedSize[0]-highPad[0] &&
		y>=enhancedSize[1]-highPad[1] && y<enhancedSize[1] &&
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
		return 8;
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
	if( x>=0 && x<lowPad[0] &&
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
		z>=enhancedSize[2]-highPad[2] && z<enhancedSize[2] )
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
// makes DIY datatype for sending and receiving one item
//
// dtype: pointer to the datatype
//
void
ItemDtype( DIY_Datatype *dtype ) {

  struct map_block_t map[1] =
  {
    {DIY_INT, OFST, 1, 0},
  };

  DIY_Create_struct_datatype(0, 1, map, dtype);

}

//
// makes MPI datatype for receiving one item
//
// cts: pointer to counts message
//
// side effects: allocates MPI datatype
//
// returns: pointer to MPI datatype
//
/*
DIY_Datatype*
RecvItemType( int *cts )
{
	DIY_Datatype *dtype = (DIY_Datatype *) malloc( sizeof( DIY_Datatype ) );

	struct map_block_t map[1] =
	{
		{DIY_FLOAT, OFST, (maxGhostSize + 4), 0 },
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
DIY_Datatype*
SendItemType( int *cts, char** pts )
{
	DIY_Datatype *dtype = (DIY_Datatype *) malloc( sizeof( DIY_Datatype ) );

	struct map_block_t map[1] =
	{
		{DIY_FLOAT, OFST, (maxGhostSize + 4), 0},
	};

	DIY_Create_struct_datatype(0, 1, map, dtype);

	return dtype;
}
*/

void*
CreateWriteType( void *item, int lid, DIY_Datatype *dtype )
{
	int min[3], size[3]; // block extents
	//DIY_Block_starts_sizes( lid, min, size );
	//int block_size = size[0] * size[1] * size[2];

	#ifdef DEBUG_MODE
	fprintf( stderr, "Block id: %d and size: %d\n", lid, block_size );
	#endif

	//DIY_Create_vector_datatype( block_size + 6, 1, DIY_FLOAT, dtype );

	return item;
}

