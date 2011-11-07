/**
 * @file MainIT_diy_global.cpp
 * Application program for global entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon
 * @author Teng-Yok
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
#include "ITL_globalentropy.h"

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

char* inputFieldFile = NULL; 	
char* patchFile = NULL;

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_histogram *histogram = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_globalentropy<SCALAR> *globalEntropyComputer_scalar = NULL;
ITL_globalentropy<VECTOR3> *globalEntropyComputer_vector = NULL;

int nBin = 100;

int rounds = 3; 				// Number of rounds of merging
//int kvalues[3] = {8, 8, 8}; 		// k-way merging, different options 
int kvalues[3] = {4, 4, 4};
//int kvalues[2] = {16, 1};
//int kvalues[1] = {2};
//int kvalues[1] = {2};

double execTime[5];
clock_t starttime, endtime;

//
// User-defined callback function for merging an array of items
// in this example we compute a global histogram by merging individual ones
//
// items: pointers to input items
// nitems: total number of input items
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: allocates resulting merged item
//
// returns: pointer to resulting merged item
//
char *ComputeMerge(char **items, int nitems) {

  // the result of the merge: type must match items in calling function--
  // we have casted away the compiler's type-checking and are on our own
  int *res = new int[nBin];

  memset(res, 0, nBin * sizeof(int));
  for (int i = 0; i < nitems; i++) {
    for (int j = 0; j < nBin; j++)
      res[j] += ((int **)items)[i][j];
  }

  return(char*)res; // cast back to char* before returning

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
char *CreateItem(int *hdr) {

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
void DeleteItem(char *item) {

  delete[] item;

}
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
void *CreateType(void *item, MPI_Datatype *dtype) {

  //MPI_Datatype *dtype = new MPI_Datatype; // DIY will free this resource for you
  //MPI_Type_contiguous(nBin, MPI_INT, dtype);
  //MPI_Type_commit(dtype); // DIY will free this resource for you
  //abs_addr = false;
  //return dtype;
  struct map_block_t map[1] = {
	{MPI_INT, OFST, nBin, 0, 1},
   };	
  DIY_Create_datatype(DIY_Addr(item), 1, map, dtype);
  
  return MPI_BOTTOM;
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

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Create the blocking and default assignment
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int given[3] = {0, 0, 0};

	// Decompose domain
	// The blocks do not need to share face in this case 
	DIY_Decompose( 0, 0, 0, given );

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
		DIY_Block_starts_sizes(i, &diy_min[3*i], &diy_size[3*i] );
	
		// post a read for the block
		if( fieldType == 0 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(data[i]));
		if( fieldType == 1 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, complex, (void**)&(vectordata[i]));
		//if( fieldType == 1 ) DIY_Add_block_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, MPI_FLOAT, (void**)&(vectordata[i]));		
		
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
	DIY_Read_all();
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );

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
			scalarField = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );
	
			// Initialize class that can compute entropy
			globalEntropyComputer_scalar = new ITL_globalentropy<SCALAR>( scalarField, histogram );

			if( histogramLowEnd != histogramHighEnd )
				globalEntropyComputer_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );				

			// Create bin field (Local analysis)
			globalEntropyComputer_scalar->computeHistogramBinField( "scalar", nBin );

			// Compute frequencies
			globalEntropyComputer_scalar->computeHistogramFrequencies( nBin );

			// Get histogram frequencies
			freqList[k] = new int[nBin];
			globalEntropyComputer_scalar->getHistogramFrequencies( nBin, freqList[k] );
		        //for (int i = 0; i < nBin; i++)
	     		//	fprintf(stderr, "freq[%d] = %d\n", i, freqList[k][i]);

			// Clear up
			delete globalEntropyComputer_scalar;
			delete scalarField;
		}
		if( fieldType == 1 )
		{
		
			// Initialize ITL vector field with block data
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

			// Initialize class that can compute entropy
			globalEntropyComputer_vector = new ITL_globalentropy<VECTOR3>( vectorField, histogram );

			// Create bin field (local analysis)
			globalEntropyComputer_vector->computeHistogramBinField( "vector", nBin );

			// Compute frequencies
			globalEntropyComputer_vector->computeHistogramFrequencies( nBin );

			// Get histogram frequencies
			freqList[k] = new int[nBin];	
			globalEntropyComputer_vector->getHistogramFrequencies( nBin, freqList[k] );

			// Clear up
			delete globalEntropyComputer_vector;
			delete vectorField; 
		}

	}// End for loop

	// Merge the local analyses
	DIY_Merge_blocks( (char**)freqList, (int **)NULL, nblocks, (char***)&hist, rounds, kvalues, &ComputeMerge, &CreateItem, &DeleteItem, &CreateType, &nb_merged );

	//for (int b = 0; b < nb_merged; b++) {
	// for (int i = 0; i < nBin; i++)
     	//	fprintf(stderr, "hist[%d][%d] = %d\n", b, i, hist[b][i]);
  	//}

	// Compute global entropy of the merged blocks
	if (nb_merged)
	{		
		for( int k=0; k<nb_merged; k++ )
		{			

			if( verboseMode == 1 )
			{			
				for( int p=0; p<nBin;p++ )
					printf( "%d, ", hist[k][p] );
				printf( "\n" );
			}
	
			// Directly compute entropy
			globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( hist[k], dataSize[0]*dataSize[1]*dataSize[2], nBin, false );

			// Print entropy	
			//printf( "Global Entropy: %f\n", globalEntropy );
		}
	}
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation Time: %f, %f seconds\n", rank, execTime[0], execTime[1] );
	else 			printf( "%d, %f, %f\n", rank, execTime[0], execTime[1] );
	
	
	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	for( int k=0; k<nblocks; k++ )
		delete [] freqList[k];
	delete [] freqList;

	// Finalize BIL
	if( fieldType == 1 ) MPI_Type_free( &complex );
	BIL_Finalize();

	// Finalize MPI
	MPI_Finalize();

}// End main




