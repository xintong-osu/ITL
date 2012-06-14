/**
 * @file MainIT_diy_crossvalidate.cpp
 * Application program for cross-validation
 * Created on: April 12, 2012
 * @author Abon
 * @author Teng-Yok
 */
#define DISPLAY_RESULTS
//#define BYTE_SWAP
#define DATA_TO_HIST_GEODESICMAPPING

#include <mpi.h>

#ifdef DISPLAY_RESULTS
#include "engine.h"
#endif

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

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

int numProcs;
int rank;

int functionType = 0;
int fieldType = 0;
int nDim = 3;
int nDivision = 4;
int nBinOptimal = 100;
int nBinStart = 1000;
int nBinCurrent = 1000;
int nIter = 100;
int step = 10;
SCALAR histogramLowEnd = 0, histogramHighEnd = 0;

int dataSize[3];
int low[3];
int high[3];
float lowF[] = {0, 0, 0};
float highF[] = {0, 0, 0};
int tot_blocks = 64;

int verboseMode = 0;

char* inputFieldFile = NULL;
char* outFile = NULL;
char* patchFile = NULL;

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_field_regular<int> *binField = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR>* histMapper_scalar = NULL;
ITL_histogrammapper<VECTOR3>* histMapper_vector = NULL;

int rounds = 3; 				// Number of rounds of merging
int kvalues[5] = {8, 8, 8, 8, 8}; 		// k-way merging, different options
//int kvalues[3] = {4, 4, 4};
//int kvalues[2] = {16, 1};
//int kvalues[1] = {2};
//int kvalues[1] = {2};

double execTime[5];
clock_t starttime, endtime;

#ifdef DISPLAY_RESULTS
Engine *engine = NULL;
mxArray *mxNBinList = NULL, *mxScoreList = NULL, *mxOptNBinList = NULL;
#endif

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
void
ComputeMerge(char **items, int* gids, int nitems )
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
		for (int j = 0; j < nBinCurrent; j++)
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

  int *bins = new int[nBinCurrent]; // DIY will free this resource for you
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
DeleteItem(char *item)
{
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
void*
CreateType(void *item, DIY_Datatype *dtype )
{

	//MPI_Datatype *dtype = new MPI_Datatype; // DIY will free this resource for you
	//MPI_Type_contiguous(nBin, MPI_INT, dtype);
	//MPI_Type_commit(dtype); // DIY will free this resource for you
	//abs_addr = false;
	//return dtype;
	//struct map_block_t map[1] = {
	//	{MPI_INT, OFST, nBin, 0, 1},
	// };
  	//DIY_Create_datatype(DIY_Addr(item), 1, map, dtype);
  	//return MPI_BOTTOM;

	DIY_Create_vector_datatype( nBinCurrent, 1, DIY_INT, dtype );
	return item;
}

/**
 * Serial cross validation.
 */
void crossvalidate_serial();
/**
 * Serial cross validation with  single data access.
 */
void crossvalidate_serial_optimized();
/**
 * Parallel cross validation
 */
void crossvalidate_parallel();
/**
 * Parallel cross validation with multiple data access
 */
void crossvalidate_parallel_optimized();
/**
 *
 */
void crossvalidate_blockwise_parallel_optimized();
#ifdef DISPLAY_RESULTS
/**
 *
 */
void displayCrossValidationResult_Matlab( double* nBinList, double* scoreList );
/**
 *
 */
void displayBlockwiseCrossValidationResult_Matlab( double* optNBinList, int nblocks );
#endif

/**
 * Main function.
 * Program starts from here.
 */
int
main( int argc, char** argv )
{
	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	functionType = atoi( ITL_util<float>::getArgWithName( "functionType", &argNames, &argValues ) );
	fieldType = atoi( ITL_util<float>::getArgWithName( "fieldType", &argNames, &argValues ) );
	inputFieldFile = ITL_util<float>::getArgWithName( "inputField", &argNames, &argValues );
	//
	nDim = atoi(  ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	dataSize[0] = atoi( ITL_util<float>::getArgWithName( "nX", &argNames, &argValues ) );
	dataSize[1] = atoi( ITL_util<float>::getArgWithName( "nY", &argNames, &argValues ) );
	dataSize[2] = atoi( ITL_util<float>::getArgWithName( "nZ", &argNames, &argValues ) );
	//
	nBinStart = atoi(  ITL_util<float>::getArgWithName( "nBinStart", &argNames, &argValues ) );
	nDivision = atoi(  ITL_util<float>::getArgWithName( "nDivision", &argNames, &argValues ) );
	nIter = atoi(  ITL_util<float>::getArgWithName( "nIter", &argNames, &argValues ) );
	step = atoi(  ITL_util<float>::getArgWithName( "step", &argNames, &argValues ) );
	//
	histogramLowEnd = (float)atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramHighEnd = (float)atof( ITL_util<float>::getArgWithName( "histHigh", &argNames, &argValues ) );
	//
	tot_blocks = atoi( ITL_util<float>::getArgWithName( "nBlock", &argNames, &argValues ) );
	//
	rounds = atoi( ITL_util<float>::getArgWithName( "numround", &argNames, &argValues ) );
	kvalues[0] = atoi( ITL_util<float>::getArgWithName( "r0", &argNames, &argValues ) );
	kvalues[1] = atoi( ITL_util<float>::getArgWithName( "r1", &argNames, &argValues ) );
	kvalues[2] = atoi( ITL_util<float>::getArgWithName( "r2", &argNames, &argValues ) );
	kvalues[3] = atoi( ITL_util<float>::getArgWithName( "r3", &argNames, &argValues ) );
	kvalues[4] = atoi( ITL_util<float>::getArgWithName( "r4", &argNames, &argValues ) );
	//
	verboseMode = atoi( ITL_util<float>::getArgWithName( "verbose", &argNames, &argValues ) );

	// Initialize ITL
	ITL_base::ITL_init();

	// Initialize histogram
	histogram = new ITL_histogram( "!" );

	// Initialize data-to-histogram converter class
	if( fieldType == 0 )
		histMapper_scalar = new ITL_histogrammapper<SCALAR>( histogram );
	else if( fieldType == 1 )
		histMapper_vector = new ITL_histogrammapper<VECTOR3>( histogram );

	switch( functionType )
	{
	case 0:
		printf( "Entering serial cross-validation ...\n" );
		crossvalidate_serial();
		break;
	case 1:
		printf( "Entering serial cross-validation with single data access ...\n" );
		crossvalidate_serial_optimized();
		break;
	case 2:
		printf( "Entering parallel cross-validation ...\n" );
		crossvalidate_parallel();
		break;
	case 3:
		printf( "Entering parallel cross-validation with single data access ...\n" );
		crossvalidate_parallel_optimized();
		break;
	case 4:
		printf( "Entering parallel cross-validation for each block ...\n" );
		crossvalidate_blockwise_parallel_optimized();
		break;

	default:
		break;
	}// end switch

	// Finalize MPI
	MPI_Finalize();
}

/**
 * Serial cross validation.
 */
void
crossvalidate_serial()
{
	int *nBinArray = new int[nIter];
	double *scoreList = new double[nIter];

	#ifdef BYTE_SWAP
	int nDataElem = dataSize[0]*dataSize[1]*dataSize[2];
	if( fieldType == 0 )
	{
		swap((char *)scalarFieldData, nDataElem, sizeof( SCALAR ));
	}
	else
	{
		swap((char *)vectorFieldData, nDataElem*3, sizeof( SCALAR ));
	}
	#endif

	starttime = ITL_util<float>::startTimer();
	if( fieldType == 0 )
	{
		// Read field
		if( verboseMode == 1 ) printf( "Reading scalar field ...\n" );
		scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinarySerial( inputFieldFile, nDim, dataSize );
		if( verboseMode == 1 ) printf( "Field dimension: %d %d %d\n", dataSize[0], dataSize[1], dataSize[2] );

		// Create a scalar field class from the file
		highF[0] = dataSize[0]-1.0f;
		highF[1] = dataSize[1]-1.0f;
		highF[2] = dataSize[2]-1.0f;
		scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF );
		if( verboseMode == 1 ) printf( "Scalar field created\n" );

		// Run cross validation
		if( verboseMode == 1 ) printf( "Starting cross-validation with %d iterations, from %d bins with step %d ...\n", nIter, nBinStart, step );
		if( histogramLowEnd != histogramHighEnd )
			histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
		histMapper_scalar->crossValidate_scalar( scalarField,
										  	  	 nIter, nBinStart, step,
										  nBinArray, scoreList,
										  &nBinOptimal );
	}
	else if( fieldType == 1 )
	{
		printf( "This kind of cross-validation is not implemented for vector field yet ... \n" );
		exit(0);
	}
	execTime[0] = ITL_util<float>::endTimer( starttime );

	printf( "Optimal number of bins: %d\n", nBinOptimal );
	printf( "Cross-validation time: %g seconds\n", execTime[0] );

	// Display
	#ifdef DISPLAY_RESULTS
	double* nBinArray_double = new double[nIter];
	double* scoreArray_double = new double[nIter];
	for( int i=0; i<nIter; i++ )
	{
		nBinArray_double[i] = nBinArray[i];
		scoreArray_double[i] = scoreList[i];
	}
	displayCrossValidationResult_Matlab( nBinArray_double, scoreArray_double );
	delete [] nBinArray_double;
	delete [] scoreArray_double;
	#endif

	// Clear up
	delete [] nBinArray;
	delete [] scoreList;
}
/**
 * Serial cross validation with  single data access.
 */
void
crossvalidate_serial_optimized()
{
	int *nBinArray = NULL;
	int *freqList = NULL;
	double *scoreList = NULL;

	#ifdef BYTE_SWAP
	int nDataElem = dataSize[0]*dataSize[1]*dataSize[2];
	if( fieldType == 0 )
	{
		swap((char *)scalarFieldData, nDataElem, sizeof( SCALAR ));
	}
	else
	{
		swap((char *)vectorFieldData, nDataElem*3, sizeof( SCALAR ));
	}
	#endif

	starttime = ITL_util<float>::startTimer();
	if( fieldType == 0 )
	{
		nBinArray = new int[nIter];
		scoreList = new double[nIter];

		// Read field
		if( verboseMode == 1 ) printf( "Reading scalar field ...\n" );
		scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinarySerial( inputFieldFile, nDim, dataSize );
		if( verboseMode == 1 ) printf( "Field dimension: %d %d %d\n", dataSize[0], dataSize[1], dataSize[2] );

		// Create a scalar field class from the file
		highF[0] = dataSize[0]-1.0f;
		highF[1] = dataSize[1]-1.0f;
		highF[2] = dataSize[2]-1.0f;
		scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF );
		if( verboseMode == 1 ) printf( "Scalar field created\n" );

		// Run cross validation
		if( verboseMode == 1 ) printf( "Starting optimized cross-validation with %d iterations, from %d bins ...\n", nIter, nBinStart );
		if( histogramLowEnd != histogramHighEnd )
			histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
		histMapper_scalar->crossValidate_optimized_scalar( scalarField,
														   nIter, nBinStart,
														   nBinArray, scoreList,
														   &nBinOptimal );
	}
	else
	{
		nBinArray = new int[nDivision];
		scoreList = new double[nDivision];

		// Read field
		if( verboseMode == 1 ) printf( "Reading vector field ...\n" );
		vectorFieldData = ITL_ioutil<VECTOR3>::readFieldBinarySerial( inputFieldFile, nDim, dataSize );
		if( verboseMode == 1 ) printf( "Field dimension: %d %d %d\n", dataSize[0], dataSize[1], dataSize[2] );

		// Create a scalar field class from the file
		highF[0] = dataSize[0]-1.0f;
		highF[1] = dataSize[1]-1.0f;
		highF[2] = dataSize[2]-1.0f;
		vectorField = new ITL_field_regular<VECTOR3>( vectorFieldData, nDim, lowF, highF );

		// Run cross validation
		if( verboseMode == 1 ) printf( "Starting optimized cross-validation with %d-level geodesic grid ...\n", nDivision );
		#ifdef DATA_TO_HIST_GEODESICMAPPING
		// nDivision: 1 -> nBin: 20
		// nDivision: 2 -> nBin: 80
		// nDivision: 3 -> nBin: 320
		// nDivision: 4 -> nBin: 1280
		histMapper_vector->createSphericalGrid( nDivision, &nBinStart );
		histMapper_vector->computeTable1( nBinStart, nDivision );
		histMapper_vector->computeHistogramBinField_Vector_Geodesic( vectorField, &binField, nBinStart, nDivision );
		freqList = new int[nBinStart];
		histMapper_vector->computeHistogramFrequencies( &binField, freqList, nBinStart );
		//for( int l=0; l<nBinStart; l++ )
		//	printf( "%d, ", freqList[l] );
		//printf( "\n" );
		#else
		printf( "Geodesic mapping not enabled ... exiting!\n" ).
		exit(0);
		#endif
		histMapper_vector->crossValidate_optimized_vector( freqList,
														   dataSize[0]*dataSize[1]*dataSize[2],
														   nBinStart, nDivision,
														   nBinArray, scoreList,
														   &nBinOptimal );
	}
	execTime[0] = ITL_util<float>::endTimer( starttime );

	printf( "Optimal number of bins: %d\n", nBinOptimal );
	printf( "Cross-validation time: %g seconds\n", execTime[0] );

	// Display
	#ifdef DISPLAY_RESULTS
	if( fieldType == 1 )
		nIter = nDivision;
	double* nBinArray_double = new double[nIter];
	double* scoreArray_double = new double[nIter];
	for( int i=0; i<nIter; i++ )
	{
		//printf( "%d %g\n", nBinArray[i], scoreList[i] );
		nBinArray_double[i] = nBinArray[i];//nIter-1-i];
		scoreArray_double[i] = scoreList[i];//[nIter-1-i];
	}
	displayCrossValidationResult_Matlab( nBinArray_double, scoreArray_double );
	delete [] nBinArray_double;
	delete [] scoreArray_double;
	#endif

	// Clear up
	delete [] nBinArray;
	delete [] scoreList;
	delete [] freqList;

}
/**
 * Parallel cross validation
 */
void
crossvalidate_parallel()
{
	int nblocks;						// My local number of blocks
	int nb_merged; 						// Number of output merged blocks

	starttime = ITL_util<float>::startTimer();

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Create the blocking and default assignment
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int given[3] = {0, 0, 0};

	// Decompose domain
	DIY_Decompose( 0, 0, 0, given );

	// Create data type for vectors
	MPI_Datatype complex;
	if( fieldType == 1 )
	{
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  		MPI_Type_commit( &complex );
	}

	// Allocate memory for pointers that will hold block data
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
	DIY_Read_all();
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );

	starttime = ITL_util<float>::startTimer();

	#ifdef BYTE_SWAP
	int nBlockElem = 0;
	for( int k=0; k<nblocks; k++ )
	{
		nBlockElem = diy_size[3*k]*diy_size[3*k+1]*diy_size[3*k+2];
		if( fieldType == 0 )
			swap((char *)&data[k][0], nBlockElem, sizeof( SCALAR ));
		else if( fieldType == 1 )
			swap((char *)&vectordata[k][0], nBlockElem*3, sizeof( SCALAR ));
	}
	#endif


  	int N = dataSize[0]*dataSize[1]*dataSize[2];
	double h;
  	double *normFreqListPtr = NULL;
	int *nBinArray = new int[nIter];
	double *scoreList = new double[nIter];
	double minScore = 10000000;
	int minScoreIndex = -1;
 	int **freqList = NULL; 	// global histogram(s), in this example (can be > 1, depending
	            			// on the degree of merging, hence the double pointer)

	// Set global range of histogram
	if( histogramLowEnd != histogramHighEnd )
		histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );

	// Compute number of bins for different iterations
	for( int i = 0; i<nIter; i++ )
		nBinArray[i] = nBinStart + i*step;

	// Iterate through different number of bins
	for( int i = 0; i<nIter; i++ )
	{
		// Allocate memory for storing frequencies
	  	freqList = new int*[nblocks];

		// Iterate over each block and compute histogram
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

				// Create bin field (Local analysis)
				//cout << "-2" << endl;
				histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBinArray[i] );

				// Compute frequencies
				//cout << "-1" << endl;
				freqList[k] = new int[nBinArray[i]];
				histMapper_scalar->computeHistogramFrequencies( &binField, freqList[k], nBinArray[i] );

				//for( int p=0; p<nBinArray[i]; p++ )
				//	printf( "%d, ", freqList[k][p] );
				//printf( "\n" );

				// Clear up
				//cout << "0" << endl;
				delete binField;
				binField = NULL;
				delete scalarField;

			}
			if( fieldType == 1 )
			{
				printf( "This kind of cross-validation is not implemented for vector field yet ... \n" );
				exit(0);
			}

		}// End for loop

		// Merge the local analyses
		nBinCurrent = nBinArray[i]; 					/*********************/
		DIY_Merge_blocks( (char**)freqList, (int **)NULL,
						  rounds, kvalues,
						  &ComputeMerge,
						  &CreateItem,
						  &CreateType, &nb_merged );

		// Compute cross-validation score based on the merged frequencies
		if (nb_merged)
		{
			//cout << "adding frequencies " << nb_merged << endl;
			for( int k=0; k<nb_merged; k++ )
			{
				//if( verboseMode == 1 )
				//{
				//	for( int p=0; p<nBinArray[i]; p++ )
				//		printf( "%d, ", freqList[k][p] );
				//	printf( "\n" );
				//}

				// Normalize merged frequencies
				normFreqListPtr = new double[nBinArray[i]];
				for( int j=0; j<nBinArray[i]; j++ )
					normFreqListPtr[j] =  freqList[k][j] / (double)N;

				// Update bin width
				h = ( histogramHighEnd - histogramLowEnd ) / (double) nBinArray[i];

				// Compute cross-validation score
				if( fieldType == 0 )
					scoreList[i] = histMapper_scalar->computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );
				else
					scoreList[i] = histMapper_vector->computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

				//#ifdef DEBUG_MODE
				printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
				//#endif

				if( scoreList[i] < minScore )
				{
					minScore = scoreList[i];
					minScoreIndex = i;
				}

				//cout << "-2" << endl;
				delete [] normFreqListPtr;

			}// end for k : merged blocks

		}// end if : nb_merged

		for( int k=0; k<nblocks; k++ )
			delete [] freqList[k];
		delete [] freqList;

	}// end for i : number of different bins

	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Print Optimal number of bins
	if( nb_merged )
	{
		// Display
		// Display
		#ifdef DISPLAY_RESULTS
		double* nBinArray_double = new double[nIter];
		double* scoreArray_double = new double[nIter];
		for( int i=0; i<nIter; i++ )
		{
			nBinArray_double[i] = nBinArray[i];
			scoreArray_double[i] = scoreList[i];
		}
		displayCrossValidationResult_Matlab( nBinArray_double, scoreArray_double );
		delete [] nBinArray_double;
		delete [] scoreArray_double;
		std::cin.get();
		#endif

		nBinOptimal = nBinArray[minScoreIndex];
		if( verboseMode == 1 )
		{
			printf( "Optimal number of bins: %d\n", nBinOptimal );
			printf( "Read/Cross-validation time: %g %g seconds\n", execTime[0], execTime[1] );
		}
		else
		{
			printf( "%d, %g, %g, %d\n", rank, execTime[0], execTime[1], nBinOptimal );
		}

	}

	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	//for( int k=0; k<nblocks; k++ )
	//	delete [] freqList[k];
	//delete [] freqList;
	delete [] nBinArray;
	delete [] scoreList;

	// Finalize BIL
	if( fieldType == 1 ) MPI_Type_free( &complex );
	DIY_Finalize();
}

/**
 * Parallel cross validation with multiple data access
 */
void
crossvalidate_parallel_optimized()
{
	int nblocks;						// My local number of blocks
	int nb_merged; 						// Number of output merged blocks

	starttime = ITL_util<float>::startTimer();

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Create the blocking and default assignment
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int given[3] = {0, 0, 0};

	// Decompose domain
	DIY_Decompose( 0, 0, 0, given );

	// Create data type for vectors
	MPI_Datatype complex;
	if( fieldType == 1 )
	{
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  		MPI_Type_commit( &complex );
	}

	// Allocate memory for pointers that will hold block data
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
	DIY_Read_all();
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );

	starttime = ITL_util<float>::startTimer();

	#ifdef BYTE_SWAP
	int nBlockElem = 0;
	for( int k=0; k<nblocks; k++ )
	{
		nBlockElem = diy_size[3*k]*diy_size[3*k+1]*diy_size[3*k+2];
		if( fieldType == 0 )
			swap((char *)&data[k][0], nBlockElem, sizeof( SCALAR ));
		else if( fieldType == 1 )
			swap((char *)&vectordata[k][0], nBlockElem*3, sizeof( SCALAR ));
	}
	#endif

	// Allocate memory for storing frequencies
 	int **freqList = NULL; 	// global histogram(s), in this example (can be > 1, depending
	            		// on the degree of merging, hence the double pointer)
  	freqList = new int*[nblocks];
	int **hist;

	if( fieldType == 1 )
	{
		#ifdef DATA_TO_HIST_GEODESICMAPPING
		// nDivision: 1 -> nBin: 20
		// nDivision: 2 -> nBin: 80
		// nDivision: 3 -> nBin: 320
		// nDivision: 4 -> nBin: 1280
		histMapper_vector->createSphericalGrid( nDivision, &nBinStart );
		histMapper_vector->computeTable1( nBinStart, nDivision );
		#else
		printf( "Geodesic mapping not enabled ... exiting!\n" ).
		exit(0);
		#endif
	}

  	// Scan through the blocks, compute global entropy and list values
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

			// Create bin field (Local analysis)
			//cout << "-2" << endl;
			if( histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBinStart );

			// Compute frequencies
			//cout << "-1" << endl;
			freqList[k] = new int[nBinStart];
			histMapper_scalar->computeHistogramFrequencies( &binField, freqList[k], nBinStart );

			// Clear up
			//cout << "0" << endl;
			delete binField;
			binField = NULL;
			delete scalarField;

		}
		if( fieldType == 1 )
		{
			// Initialize ITL vector field with block data
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

			// Create bin field (local analysis)
			if( verboseMode == 1 ) printf( "Starting parallel optimized cross-validation with %d-level geodesic grid ...\n", nDivision );
			#ifdef DATA_TO_HIST_GEODESICMAPPING
			histMapper_vector->computeHistogramBinField_Vector_Geodesic( vectorField, &binField, nBinStart, nDivision );
			#else
			printf( "Geodesic mapping not enabled ... exiting!\n" ).
			exit(0);
			#endif

			// Get histogram frequencies
			freqList[k] = new int[nBinStart];
			histMapper_vector->computeHistogramFrequencies( &binField, freqList[k], nBinStart );

			// Clear up
			delete binField;
			binField = NULL;
			delete vectorField;
		}

	}// End for loop

	// Merge the local analyses
	DIY_Merge_blocks( (char**)freqList, (int **)NULL,
					  rounds, kvalues,
					  &ComputeMerge,
					  &CreateItem,// &DeleteItem,
					  &CreateType, &nb_merged );
	//cout << "Merged " << endl;

	// Compute cross-validation score based on the merged frequencies
	if (nb_merged)
	{
		int *nBinArray = new int[nIter];
		double *scoreList = new double[nIter];

		assert( hist != NULL );
		//cout << "adding frequencies " << nb_merged << endl;
		for( int k=0; k<nb_merged; k++ )
		{
			//if( verboseMode == 1 )
			//{
			//	for( int p=0; p<nBinStart;p++ )
			//		printf( "%d, ", freqList[k][p] );
			//	printf( "\n" );
			//}

			// Run cross-validation using the combined frequency list
			if( fieldType == 0 )
				histMapper_scalar->crossValidate_optimized2( freqList[k], nIter, nBinStart,
															 dataSize[0]*dataSize[1]*dataSize[2],
															 nBinArray, scoreList, &nBinOptimal );
			else
				histMapper_vector->crossValidate_optimized_vector( freqList[k],
																   dataSize[0]*dataSize[1]*dataSize[2],
																   nBinStart,
																   nDivision,
																   nBinArray, scoreList,
																   &nBinOptimal );

			// Display
			#ifdef DISPLAY_RESULTS
			if( fieldType == 1 )
				nIter = nDivision;
			double* nBinArray_double = new double[nIter];
			for( int i=0; i<nIter; i++ )
				nBinArray_double[i] = nBinArray[i];
			displayCrossValidationResult_Matlab( nBinArray_double, scoreList );
			delete [] nBinArray_double;
			#endif

			execTime[1] = ITL_util<float>::endTimer( starttime );

			// Print Optimal number of bins
			if( verboseMode == 1 )
			{
				printf( "Optimal number of bins: %d\n", nBinOptimal );
				printf( "Read/Cross-validation time: %g %g seconds\n", execTime[0], execTime[1] );
			}
			else
			{
				printf( "%d, %g, %g, %d\n", rank, execTime[0], execTime[1], nBinOptimal );
			}
		}

		delete [] nBinArray;
		delete [] scoreList;
	}


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

}// end method

/**
 * Parallel cross validation with multiple data access
 */
void
crossvalidate_blockwise_parallel_optimized()
{
	int nblocks;						// My local number of blocks
	int nb_merged; 						// Number of output merged blocks

	starttime = ITL_util<float>::startTimer();

	// Initialize DIY after initializing MPI
	DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, MPI_COMM_WORLD );
	if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

	// Create the blocking and default assignment
	// Note in the blocking call that we are not adding extra ghost cells, but we
	// are sharing boundaries between blocks (share_face = 1)
	int given[3] = {0, 0, 0};

	// Decompose domain
	DIY_Decompose( 0, 0, 0, given );

	// Create data type for vectors
	MPI_Datatype complex;
	if( fieldType == 1 )
	{
	 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  		MPI_Type_commit( &complex );
	}

	// Allocate memory for pointers that will hold block data
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
	DIY_Read_all();
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 )	printf( "Data read complete .. all processes synced\n" );

	starttime = ITL_util<float>::startTimer();

	#ifdef BYTE_SWAP
	int nBlockElem = 0;
	for( int k=0; k<nblocks; k++ )
	{
		nBlockElem = diy_size[3*k]*diy_size[3*k+1]*diy_size[3*k+2];
		if( fieldType == 0 )
			swap((char *)&data[k][0], nBlockElem, sizeof( SCALAR ));
		else if( fieldType == 1 )
			swap((char *)&vectordata[k][0], nBlockElem*3, sizeof( SCALAR ));
	}
	#endif

	// Allocate memory for storing frequencies
 	int **freqList = NULL; 	// global histogram(s), in this example (can be > 1, depending
	            		// on the degree of merging, hence the double pointer)
  	freqList = new int*[nblocks];
  	int* optimalNumBinList = new int[nblocks];

  	// Scan through the blocks, compute global entropy and list values
	for( int k=0; k<nblocks; k++ )
	{
		lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
		lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
		lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

		int *nBinArray = new int[nIter];
		double *scoreList = new double[nIter];

		if( fieldType == 0 )
		{
			// Initialize ITL scalar field with block data
			//cout << "-4" << endl;
			scalarField = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );

			// Create bin field (Local analysis)
			//cout << "-2" << endl;
			if( histogramLowEnd != histogramHighEnd )
				histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
			histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBinStart );

			// Compute frequencies
			//cout << "-1" << endl;
			freqList[k] = new int[nBinStart];
			histMapper_scalar->computeHistogramFrequencies( &binField, freqList[k], nBinStart );

			// Run cross-validation
			histMapper_scalar->crossValidate_optimized2( freqList[k], nIter, nBinStart,
														 dataSize[0]*dataSize[1]*dataSize[2],
														 nBinArray, scoreList, &optimalNumBinList[k] );

			// Clear up
			//cout << "0" << endl;
			delete binField;
			binField = NULL;
			delete scalarField;

		}
		if( fieldType == 1 )
		{
			// Initialize ITL vector field with block data
			vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

			// Create bin field (local analysis)
			histMapper_vector->computeHistogramBinField_Vector( vectorField, &binField, nBinStart );

			// Get histogram frequencies
			freqList[k] = new int[nBinStart];
			histMapper_vector->computeHistogramFrequencies( &binField, freqList[k], nBinStart );

			histMapper_vector->crossValidate_optimized2( freqList[k], nIter, nBinStart,
														 dataSize[0]*dataSize[1]*dataSize[2],
														 nBinArray, scoreList, &optimalNumBinList[k] );

			// Clear up
			delete binField;
			binField = NULL;
			delete vectorField;
		}

		execTime[1] = ITL_util<float>::endTimer( starttime );

		// Print Optimal number of bins
		if( verboseMode == 1 )
		{
			printf( "%d : Optimal number of bins: %d\n", k, optimalNumBinList[k] );
			printf( "%d : Cross-validation time: %g seconds\n", k, execTime[0] );
		}
		else
		{
			printf( "%d, %d, %g, %g, %d\n", rank, k, execTime[0], execTime[1], optimalNumBinList[k] );
		}


		delete [] nBinArray;
		delete [] scoreList;

	}// End for k : blocks

	// Display
	#ifdef DISPLAY_RESULTS
	double* optNBinList_double = new double[nblocks];
	for( int i=0; i<nblocks; i++ )
		optNBinList_double[i] = optimalNumBinList[i];
	displayBlockwiseCrossValidationResult_Matlab( optNBinList_double, nblocks );
	delete [] optNBinList_double;
	#endif

	// Clear up
	delete [] diy_min;
	delete [] diy_max;
	delete [] diy_size;
	for( int k=0; k<nblocks; k++ )
		delete [] freqList[k];
	delete [] freqList;
	delete [] optimalNumBinList;

	// Finalize BIL
	if( fieldType == 1 ) MPI_Type_free( &complex );
	DIY_Finalize();

}// end method

#ifdef DISPLAY_RESULTS
void
displayCrossValidationResult_Matlab( double* nBinList, double* scoreList )
{
	// Call engOpen with a NULL string. This starts a MATLAB process
    // on the current host using the command "matlab".
	fprintf( stderr, "Initializing Matlab engine ...\n" );

	if (!(engine = engOpen("")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(1);
	}

	mxNBinList = mxCreateDoubleMatrix( 1, nIter, mxREAL );
	memcpy( (void *)mxGetPr( mxNBinList ), (void *)nBinList, sizeof(double)*(nIter) );

	mxScoreList = mxCreateDoubleMatrix( 1, nIter, mxREAL );
	memcpy( (void *)mxGetPr( mxScoreList ), (void *)scoreList, sizeof(double)*(nIter) );

	// Place the variables into the MATLAB workspace
	engPutVariable( engine, "mxNBinList", mxNBinList );
	engPutVariable( engine, "mxScoreList", mxScoreList );

	// Plot
	fprintf( stderr, "Displaying cross-validation result ...\n" );
	engEvalString( engine, "displayCrossValidationResult( mxNBinList, mxScoreList );" );

	// Pause
	 std::cin.get();

	// Close matlab engine
	engClose( engine );

}

void
displayBlockwiseCrossValidationResult_Matlab( double* optNBinList, int nblocks )
{
	// Call engOpen with a NULL string. This starts a MATLAB process
    // on the current host using the command "matlab".
	fprintf( stderr, "Initializing Matlab engine ...\n" );

	if (!(engine = engOpen("")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(1);
	}

	mxOptNBinList = mxCreateDoubleMatrix( 1, nblocks, mxREAL );
	memcpy( (void *)mxGetPr( mxOptNBinList ), (void *)optNBinList, sizeof(double)*(nblocks) );

	// Place the variables into the MATLAB workspace
	engPutVariable( engine, "mxOptNBinList", mxOptNBinList );

	// Plot
	fprintf( stderr, "Displaying blockwise cross-validation result ...\n" );
	engEvalString( engine, "displayBlockwiseCrossValidationResult( mxOptNBinList );" );

	// Pause
	 std::cin.get();

	// Close matlab engine
	engClose( engine );

}
#endif

