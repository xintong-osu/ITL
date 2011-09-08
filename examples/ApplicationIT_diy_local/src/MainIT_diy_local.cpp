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

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_localentropy<SCALAR> *localEntropyComputer_scalar = NULL;
ITL_localentropy<VECTOR3> *localEntropyComputer_vector = NULL;

int tot_blocks = 512;

double execTime[3];
clock_t starttime, endtime;

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
	MPI_Datatype complex;
	if( fieldType == 1 )
	{
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  		MPI_Type_commit( &complex );
	}

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Decompose domain (with ghost layers) 
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	DIY_Decompose( 0, 6, 0, given );

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
		if( fieldType == 0 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(data[i]));
		if( fieldType == 1 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, complex, (void**)&(vectordata[i]));

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
	if( verboseMode == 1 )	printf( "Data read complete .. al processes synced\n" );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Allocate memory for storing local entropy fields for each block 		
	float *localEntropyList[nblocks];

	// Scan through the blocks, compute local entropy fields 
	starttime = ITL_util<float>::startTimer();	
	for( int k=0; k<nblocks; k++ )
	{

		lowF[0] = 0; 		highF[0] = diy_max[3*k] - diy_min[3*k];
		lowF[1] = 0; 		highF[1] = diy_max[3*k+1] - diy_min[3*k+1];
		lowF[2] = 0; 		highF[2] = diy_max[3*k+2] - diy_min[3*k+2];

		if( fieldType == 0 )
		{
			// Initialize ITL scalar field with block data

			scalarField = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );
	
			// Initialize class that can compute entropy
			localEntropyComputer_scalar = new ITL_localentropy<SCALAR>( scalarField );

			localEntropyComputer_scalar->setHistogramRange( 0.0f, 255.0f );

			// Create bin field
			localEntropyComputer_scalar->computeHistogramBinField( "scalar", nBin );

			// Compute entropy 
			localEntropyComputer_scalar->computeEntropyOfField( nBin, false );

			// Save local entropy field
			//localEntropyList[k] = new float[scalarField->grid->nVertices + 6];
			localEntropyList[k] = new float[dataSize[0]*dataSize[1]*dataSize[2]];
			localEntropyList[k][0] = diy_min[3*k];
			localEntropyList[k][1] = diy_min[3*k+1];
			localEntropyList[k][2] = diy_min[3*k+2];
			localEntropyList[k][3] = diy_max[3*k];
			localEntropyList[k][4] = diy_max[3*k+1];
			localEntropyList[k][5] = diy_max[3*k+2];
			memcpy( localEntropyList[k]+6, localEntropyComputer_scalar->getEntropyField(), sizeof(SCALAR)*scalarField->grid->nVertices );
			if( verboseMode == 1 ) printf( "Block Limits: %d, %d, %d, %d, %d, %d, local entropy computed.\n", scalarField->grid->lowInt[0],
															scalarField->grid->highInt[0],
															scalarField->grid->lowInt[1],
															scalarField->grid->highInt[1],
															scalarField->grid->lowInt[2],
															scalarField->grid->highInt[2] );

			float max = -100000;
			float min = 100000;
			for( int i=0; i<scalarField->grid->nVertices; i++ )
			{
				if( localEntropyList[k][i+6] > max ) max = localEntropyList[k][i+6];
				if( localEntropyList[k][i+6] < min ) min = localEntropyList[k][i+6];
			}

			printf( "%f %f\n", min, max );

			// Clear up
			delete localEntropyComputer_scalar;
			delete scalarField;
		}
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

	}// End for loop
	execTime[1] = ITL_util<float>::endTimer( starttime );


	// Write local entropy field 
	if( verboseMode == 1 ) printf( "Writing local entropy field ...\n" );
	starttime = ITL_util<float>::startTimer();	
	MPI_Datatype *dtype = new MPI_Datatype; // datatype for output
	MPI_Type_contiguous( dataSize[0]*dataSize[1]*dataSize[2], MPI_FLOAT, dtype );
	MPI_Type_commit( dtype );
	MPI_Barrier( MPI_COMM_WORLD ); // everyone synchronizes again
		
	float **listptr = new float*[nblocks];
	for( int i=0; i<nblocks; i++ )
	{
		listptr[i] = localEntropyList[i];
	}

	DIY_Write_open_all( outFile, 0 );
	DIY_Write_blocks_all( (void **)listptr, nblocks, dtype );
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




