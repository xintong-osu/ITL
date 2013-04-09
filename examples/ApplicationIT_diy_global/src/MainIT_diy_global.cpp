/**
 * @file MainIT_diy_global.cpp
 * Application program for global entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon
 * @author Teng-Yok
 */

//#define BYTE_SWAP

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
#include "ITL_histogramconstants.h"
#include "ITL_histogram.h"
#include "ITL_histogrammapper.h"
#include "ITL_field_regular.h"
#include "ITL_globalentropy.h"

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

int num_threads = 4; // number of threads DIY can use

char* inputFieldFile = NULL; 	
char* patchFile = NULL;

SCALAR* scalarFieldData = NULL;
VECTOR3* vectorFieldData = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR>* histMapper_scalar = NULL;
ITL_histogrammapper<VECTOR3>* histMapper_vector = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_field_regular<int> *binField = NULL;

ITL_globalentropy<SCALAR> *globalEntropyComputer_scalar = NULL;
ITL_globalentropy<VECTOR3> *globalEntropyComputer_vector = NULL;

int nBin = 100;

int rounds = 3; 					// Number of rounds of merging
int kvalues[5] = {8, 8, 8, 8, 8}; 	// k-way merging, different options
//int kvalues[3] = {4, 4, 4};
//int kvalues[2] = {16, 1};
//int kvalues[1] = {2};
//int kvalues[1] = {2};
//int kvalues[1] = {1};

double execTime[5];
clock_t starttime, endtime;

//
// User-defined callback function for merging an array of items
// in this example we compute a global histogram by merging individual ones
//
// items: pointers to input items
// nitems: total number of input items
// char * is used as a generic pointers to bytes, not necessarily to strings
// hdr: quantity information for items[0] (unused in this example)
//
// side effects: allocates resulting merged item
//
// returns: pointer to resulting merged item
//
void
ComputeMerge(char **items, int* gids, int nitems, int* hdr)
{
	// the result of the merge: type must match items in calling function--
	// we have casted away the compiler's type-checking and are on our own
	//int *res = new int[nBin];
	//memset(res, 0, nBin * sizeof(int));
	//for (int i = 0; i < nitems; i++) {
	// for (int j = 0; j < nBin; j++)
	//    res[j] += ((int **)items)[i][j];
	//}
	//return(char*)res; // cast back to char* before returning
	for (int i = 1; i < nitems; i++)
	{
		for (int j = 0; j < nBin; j++)
			((int **)items)[0][j] += ((int **)items)[i][j];
	}
}
//
// user-defined callback function for creating a received item
//
// hdr: quantity information for allocating custom parts of the item
//  (not used in this example)
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: allocates the item
//
// returns: pointer to the item
//
char*
CreateItem(int *hdr)
{

  int *bins = new int[nBin]; // DIY will free this resource for you
  return (char *)bins;
}
//
// user-defined callback function for deleting a received item
//
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: deletes the item
//
void
DeleteItem(void *item)
{
  delete[] (int *)item;
}

/**
 * user-defined callback function for destroying a received item
 * @item item to be destroyed
 */
void DestroyItem(void *item)
{
	delete[] (int *)item;
}
//
// user-defined callback function for creating an MPI datatype for the
//   received item being merged
//
// item: pointer to the item
// dtype: pointer to the datatype
// hdr: quantity information (unused in this example)
//
// side effects: commits the DIY datatype but DIY will cleanup datatype for you
//
void
CreateMergeType(void *item, DIY_Datatype *dtype, int* hdr )
{
	DIY_Create_vector_datatype( nBin, 1, DIY_INT, dtype );
}

/**
 * Main function.
 * Program starts from here.
 */
int
main( int argc, char** argv )
{
	int numProcs;
	int rank;

	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	int nDim;
	
	int dataSize[3];

	int fieldType = -1;
	int method = 0;
	int verboseMode = 1;

	float lowF[3];
	float highF[3];
	int nPartition[3]; 
	float histogramLowEnd = 0;
	float histogramHighEnd = 0;

	int nblocks;						// My local number of blocks
	int tot_blocks = 512;					// Total number of blocks across all processors
	int nb_merged; 						// Number of output merged blocks

	float globalEntropy = 0;

	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	inputFieldFile = ITL_util<float>::getArgWithName( "inputField", &argNames, &argValues );
	patchFile =  ITL_util<float>::getArgWithName( "patchFile", &argNames, &argValues );
	fieldType = atoi( ITL_util<float>::getArgWithName( "fieldType", &argNames, &argValues ) );
	nDim = atoi( ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	nBin = atoi( ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
	rounds = atoi( ITL_util<float>::getArgWithName( "numround", &argNames, &argValues ) );
	kvalues[0] = atoi( ITL_util<float>::getArgWithName( "r0", &argNames, &argValues ) );
	kvalues[1] = atoi( ITL_util<float>::getArgWithName( "r1", &argNames, &argValues ) );
	kvalues[2] = atoi( ITL_util<float>::getArgWithName( "r2", &argNames, &argValues ) );
	kvalues[3] = atoi( ITL_util<float>::getArgWithName( "r3", &argNames, &argValues ) );
	kvalues[4] = atoi( ITL_util<float>::getArgWithName( "r4", &argNames, &argValues ) );
	dataSize[0] = atoi( ITL_util<float>::getArgWithName( "nX", &argNames, &argValues ) );
	dataSize[1] = atoi( ITL_util<float>::getArgWithName( "nY", &argNames, &argValues ) );
	dataSize[2] = atoi( ITL_util<float>::getArgWithName( "nZ", &argNames, &argValues ) );
	tot_blocks = atoi( ITL_util<float>::getArgWithName( "nBlock", &argNames, &argValues ) );
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

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, dataSize, num_threads, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Create the blocking and default assignment
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int given[3] = {0, 0, 0};

	// Decompose domain
	// The blocks do not need to share face in this case 
	int did = DIY_Decompose( ROUND_ROBIN_ORDER, tot_blocks, &nblocks, 0, 0, given );

	// Allocate memory for pointers that will hold block data
	MPI_Datatype complex; 
	if( fieldType == 1 )  	
	{		
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex ); 
  		MPI_Type_commit( &complex ); 
	}
	
	float* data[nblocks];	
	VECTOR3* vectordata[nblocks];		
	if( fieldType == 0 ) memset( data, 0, sizeof(float*) * nblocks );
	if( fieldType == 1 ) memset( vectordata, 0, sizeof(VECTOR3*) * nblocks );

	// Serially visit the blocks (?)
	int* diy_min = new int[3*nblocks];
	int* diy_max = new int[3*nblocks];
	int* diy_size = new int[3*nblocks];

	for (int i = 0; i < nblocks; i++)
	{ 
		DIY_Block_starts_sizes( did, i, &diy_min[3*i], &diy_size[3*i] );
	
		// post a read for the block
		if( fieldType == 0 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, DIY_FLOAT, (void**)&(data[i]));
		if( fieldType == 1 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, complex, (void**)&(vectordata[i]));

		// print the block bounds
		for (int j = 0; j < 3; j++)
			diy_max[3*i+j] = diy_min[3*i+j] + diy_size[3*i+j] - 1;

		if( verboseMode == 1 )
		{	printf("process rank = %d "
				"block local id = %d " 
				"min = [%d %d %d] "
				"max = [%d %d %d] "
				"size = [%d %d %d]\n", 
				rank, i, 
				diy_min[3*i], diy_min[3*i+1], diy_min[3*i+2], 
				diy_max[3*i], diy_max[3*i+1], diy_max[3*i+2], 
				diy_size[3*i], diy_size[3*i+1], diy_size[3*i+2] );
		}
	}

	// Read actual data (everyone synchronizes after reading data)
	starttime = ITL_util<float>::startTimer();	
	DIY_Read_data_all();
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );

	#ifdef BYTE_SWAP
	int nBlockElem = 0;
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

	// Allocate memory for storing frequencies
 	int **freqList = NULL; 	// global histogram(s), in this example (can be > 1, depending
	            		// on the degree of merging, hence the double pointer)
  	freqList = new int*[nblocks];
	int **hist;

  	// Scan through the blocks, compute global entropy and list values 
	starttime = ITL_util<float>::startTimer();	
	for( int k=0; k<nblocks; k++ )
	{

		lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
		lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
		lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

		if( fieldType == 0 )
		{

			// Initialize ITL scalar field with block data
			//cout << "-4" << endl;
			scalarField = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );

			//cout << "-3" << endl;
			if( histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );

			// Create bin field (Local analysis)
			//cout << "-2" << endl;
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBin );

			// Initialize class that can compute entropy
			//cout << "-1" << endl;
			globalEntropyComputer_scalar = new ITL_globalentropy<SCALAR>( binField, histogram, nBin );

			// Compute frequencies
			//cout << "0" << endl;
			globalEntropyComputer_scalar->computeHistogramFrequencies();

			// Get histogram frequencies
			//cout << "1" << endl;
			freqList[k] = new int[nBin];
			globalEntropyComputer_scalar->getHistogramFrequencies( freqList[k] );
		    //for (int i = 0; i < nBin; i++)
	     	//	fprintf(stderr, "freq[%d] = %d\n", i, freqList[k][i]);
			//cout << "2" << endl;

			// Clear up
			delete globalEntropyComputer_scalar;
			delete binField;
			binField = NULL;
			delete scalarField;

			//cout << "3" << endl;

		}
		if( fieldType == 1 )
		{
		
			// Initialize ITL vector field with block data
			//cout << "-4" << endl;
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

			// Create bin field (local analysis)
			//cout << "-3" << endl;
			histMapper_vector->computeHistogramBinField_Vector( vectorField, &binField, nBin );

			// Initialize class that can compute entropy
			//cout << "-2" << endl;
			globalEntropyComputer_vector = new ITL_globalentropy<VECTOR3>( binField, histogram, nBin );

			// Compute frequencies
			//cout << "-1" << endl;
			globalEntropyComputer_vector->computeHistogramFrequencies();

			// Get histogram frequencies
			//cout << "0" << endl;
			freqList[k] = new int[nBin];	
			globalEntropyComputer_vector->getHistogramFrequencies( freqList[k] );

			// Clear up
			delete globalEntropyComputer_vector;
			delete binField;
			binField = NULL;
			delete vectorField; 
		}

	}// End for loop

	execTime[1] = ITL_util<float>::endTimer( starttime );

	starttime = ITL_util<float>::startTimer();
	// Merge the local analyses
	DIY_Merge_blocks( did, (char**)freqList, (int **)NULL,
					  rounds, kvalues,
					  &ComputeMerge,
					  &CreateItem, &DestroyItem,
					  &CreateMergeType, &nb_merged );

	//for (int b = 0; b < nb_merged; b++) {
	// for (int i = 0; i < nBin; i++)
     	//	fprintf(stderr, "hist[%d][%d] = %d\n", b, i, hist[b][i]);
  	//}

	// Compute global entropy of the merged blocks
	if (nb_merged)
	{		
		//assert( hist != NULL );
		//cout << "adding frequencies " << nb_merged << endl;
		for( int k=0; k<nb_merged; k++ )
		{			
			//if( verboseMode == 1 )
			//{
			//	for( int p=0; p<nBin;p++ )
			//		printf( "%d, ", freqList[k][p] );
			//	printf( "\n" );
			//}
	
			// Directly compute entropy
			//cout << dataSize[0] << " " << dataSize[1] << " " << dataSize[2] << endl;
			//globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( freqList[k], dataSize[0]*dataSize[1]*dataSize[2], nBin, false );

			// Print entropy	
			printf( "Global Entropy: %f\n", globalEntropy );
		}
	}
	execTime[2] = ITL_util<float>::endTimer( starttime );

	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation/Communication Time: %f, %f seconds\n", rank, execTime[0], execTime[1], execTime[2] );
	else					printf( "%d, %g, %g, %g\n", rank, execTime[0], execTime[1], execTime[2] );

	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	for( int k=0; k<nblocks; k++ )
		delete [] freqList[k];
	delete [] freqList;

	// Finalize BIL
	if( fieldType == 1 ) MPI_Type_free( &complex );
	DIY_Finalize();

	// Finalize MPI
	MPI_Finalize();

}// End main




