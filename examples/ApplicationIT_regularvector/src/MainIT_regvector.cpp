/**
 * @file MainIT_regvector.cpp
 * Application program for entropy computation of a vector field.
 * Created on: Nov 17, 2010
 * @author Abon
 * @author Teng-Yok
 */

#include <mpi.h>

#include "ITL_header.h"
#include "ITL_base.h"
#include "ITL_util.h"
#include "ITL_ioutil.h"
#include "ITL_vectormatrix.h"
#include "ITL_histogramconstants.h"
#include "ITL_localentropy.h"
#include "ITL_globalentropy.h"
#include "ITL_localjointentropy.h"
#include "ITL_field_regular.h"

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

int functionType = 0;
int sizeNeighborhood = 2;
int sizeNeighborhoodArray[3];
int nBin = 360;
int nDim = 3;
int verboseMode = 0;
int numProcs;
int myId;

int dataDim[3];
int blockDim[3];
int nBlock[3];
int blockId[3];
int low[3];
int high[3];
float lowF[] = {0, 0, 0};
float highF[] = {0, 0, 0};
int lowPad[] = {0, 0, 0};
int highPad[] = {0, 0, 0};
double execTime[5];
clock_t starttime, endtime;

VECTOR3* data = NULL;
VECTOR3* data2 = NULL;

const char *vectorFieldFile = NULL;
const char *vectorFieldFile2 = NULL;
const char *patchFile = NULL;
const char *outFieldFile = NULL;

ITL_histogram *histogram = NULL;

ITL_field_regular<VECTOR3> *vectorField = NULL;
ITL_field_regular<VECTOR3> *vectorField2 = NULL;

ITL_localentropy<VECTOR3> *localEntropyComputer = NULL;
ITL_globalentropy<VECTOR3> *globalEntropyComputer = NULL;
ITL_localjointentropy<VECTOR3> *jointEntropyComputer = NULL;

/**
 * Serial global entropy computation function with runtime computation.
 *
 */
void compute_globalentropy_serial();
/**
 * Parallel global entropy computation function for regular vector field.
 */
void compute_globalentropy_parallel();
/**
 * Parallel blockwise global entropy computation function for regular vector field.
 */
void compute_blockwiseglobalentropy_parallel();
/**
 * Serial entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_localentropy_serial();
/**
 * Parallel entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_localentropy_parallel();
/**
 * Serial joint entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_jointlocalentropy_serial();
/**
 * Parallel joint entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_jointlocalentropy_parallel();

/**
 * Main function.
 * Program starts from here.
 */
int main( int argc, char** argv )
{
	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myId );
	
	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	functionType = atoi( ITL_util<float>::getArgWithName( "functionType", &argNames, &argValues ) );
	vectorFieldFile = ITL_util<float>::getArgWithName( "vectorField", &argNames, &argValues );
	vectorFieldFile2 = ITL_util<float>::getArgWithName( "vectorField2", &argNames, &argValues );
	patchFile = ITL_util<float>::getArgWithName( "patchFile", &argNames, &argValues );
	outFieldFile = 	ITL_util<float>::getArgWithName( "outField", &argNames, &argValues );
	nDim = atoi( ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues  ) );
	nBin = atoi( ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
	sizeNeighborhood = atoi( ITL_util<float>::getArgWithName( "neighborhoodSize", &argNames, &argValues ) );
	sizeNeighborhoodArray[0] = atoi( ITL_util<float>::getArgWithName( "neighborhoodSizeX", &argNames, &argValues ) );
	sizeNeighborhoodArray[1] = atoi( ITL_util<float>::getArgWithName( "neighborhoodSizeY", &argNames, &argValues ) );  
	sizeNeighborhoodArray[2] = atoi( ITL_util<float>::getArgWithName( "neighborhoodSizeZ", &argNames, &argValues ) );  
	nBlock[0] = atoi( ITL_util<float>::getArgWithName( "nBlockX", &argNames, &argValues ) );
	nBlock[1] = atoi( ITL_util<float>::getArgWithName( "nBlockY", &argNames, &argValues ) );
	nBlock[2] = atoi( ITL_util<float>::getArgWithName( "nBlockZ", &argNames, &argValues ) );
	verboseMode = atoi( ITL_util<float>::getArgWithName( "verbose", &argNames, &argValues ) );
	
	// Initialize ITL
	ITL_base::ITL_init();

	// Initialize histogram
	histogram = new ITL_histogram( patchFile, nBin );

	switch( functionType )
	{
	case 0:
		printf( "Entering serial computation of global entropy field from a regular vector field ...\n" );
		compute_globalentropy_serial();
		break;
	case 1:
		printf( "Entering parallel computation of global entropy field from a regular vector field ...\n" );
		compute_globalentropy_parallel();
		break;
	case 2:
		printf( "Entering parallel computation of blockwise global entropy field from a regular vector field ...\n" );
		compute_blockwiseglobalentropy_parallel();
		break;
	case 3:
		printf( "Entering serial computation of local entropy field from the regular vector field ...\n" );	
		compute_localentropy_serial();
		break;
	case 4:
		printf( "Entering parallel computation of local entropy field from the regular vector field ...\n" );	
		compute_localentropy_parallel();
		break;
	case 5:
		printf( "Entering serial computation of local joint entropy field from two regular vector fields ...\n" );	
		compute_jointlocalentropy_serial();
		break;
	case 6:
		printf( "Entering parallel computation of local joint entropy fields from two regular vector fields ...\n" );	
		compute_jointlocalentropy_parallel();
		break;
	default:
		printf( "Error: Unknown Input parameter...\n" );
		break;
	}// end switch

	// Clean up
	if( localEntropyComputer != NULL ) delete localEntropyComputer;
	if( globalEntropyComputer != NULL ) delete globalEntropyComputer;
	if( jointEntropyComputer != NULL ) delete jointEntropyComputer;
	if( vectorField != NULL ) delete vectorField;
	if( vectorField2 != NULL ) delete vectorField2;

	// Finalize MPI
	MPI_Finalize();
	
}// end main

void create_raw_file( int nArg, char** argV )
{
	ITL_ioutil<VECTOR3>::truncateFileHeader( argV[3], argV[4], nDim );
}// end function

void compute_globalentropy_serial()
{
	// Read chunk of data from file
	if( verboseMode == 1 ) printf( "Reading vector field ...\n" );		
	starttime = ITL_util<float>::startTimer();	
	data = ITL_ioutil<VECTOR3>::readFieldBinarySerial( vectorFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read vector field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );

	// Create a vector field class from the file
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim, lowF, highF );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<VECTOR3>( vectorField, histogram );

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins ...\n" );
	starttime = ITL_util<float>::startTimer();
	globalEntropyComputer->computeHistogramBinField( "vector", nBin );

	if( verboseMode == 1 )
	{
		int freqList[nBin];
		globalEntropyComputer->computeHistogramFrequencies( nBin );
		globalEntropyComputer->getHistogramFrequencies( nBin, freqList );
		for( int i=0; i<nBin; i++ )
			printf( "f[%d] = %d\t", i, freqList[i] );
	}

	// Global entropy Computation
	if( verboseMode == 1 ) printf( "Computing entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	globalEntropyComputer->computeGlobalEntropyOfField( nBin, false );
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Print global entropy
	printf( "Global entropy of vector field: %f\n", globalEntropyComputer->getGlobalEntropy() );

	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation Time: %f, %f seconds\n", myId, execTime[0], execTime[1] );
	else 			printf( "%d, %f, %f\n", myId, execTime[0], execTime[1] );

}// end function

void compute_globalentropy_parallel()
{
	// Read chunk of data from file
	if( verboseMode == 1 ) printf( "Reading vector field ...\n" );	
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel2( vectorFieldFile, nDim, dataDim,
								blockDim, nBlock,blockId,
								low, high,
								lowPad, highPad,
								sizeNeighborhood, myId, numProcs);
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read vector field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );

	// Create a vector field class from the file
	highF[0] = blockDim[0]-1.0f;
	highF[1] = blockDim[1]-1.0f;
	highF[2] = blockDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim, lowF, highF );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<VECTOR3>( vectorField, histogram );

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins at each point of the vector field ...\n" ); 
	starttime = ITL_util<float>::startTimer();
	globalEntropyComputer->computeHistogramBinField( "vector", nBin );

	// Compute frequencies
	globalEntropyComputer->computeHistogramFrequencies( nBin );

	// Get histogram frequencies
	#if defined( _WIN32 ) || defined( _WIN64 )
		int* freqList = new int[nBin];
	#else
		int freqList[nBin];
	#endif
	globalEntropyComputer->getHistogramFrequencies( nBin, freqList );
		
	// Sync with all processors and and sum up all histogram frequencies to processor 0
	MPI_Barrier( MPI_COMM_WORLD );

	#if defined( _WIN32 ) || defined( _WIN64 )
		int* reducedFreqList = new int[nBin];
	#else
		int reducedFreqList[nBin];
	#endif
	MPI_Reduce( freqList, reducedFreqList, nBin, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );  

	// Compute and print global entropy	
	if( myId == 0 )	
	{
		float globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( reducedFreqList, dataDim[0]*dataDim[1]*dataDim[2], nBin, false );
		printf( "Global entropy of vector field: %f\n", globalEntropy );
		if( verboseMode == 1 )
		{
			for( int i=0; i<nBin; i++ )
				printf( "%d:f_reduced[%d]=%d\n", myId, i, reducedFreqList[i] );
		}	
	}
	execTime[1] = ITL_util<float>::endTimer( starttime );
	
	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation Time: %f, %f seconds\n", myId, execTime[0], execTime[1] );
	else 			printf( "%d, %f, %f\n", myId, execTime[0], execTime[1] );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] freqList;
		delete [] reducedFreqList;
	#endif
	
}// end function

void compute_blockwiseglobalentropy_parallel()
{
	// Read chunk of data from file
	if( verboseMode == 1 ) printf( "Reading vector field ...\n" );	
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel2( vectorFieldFile, nDim, dataDim,
								blockDim, nBlock,blockId,
								low, high,
								lowPad, highPad,
								sizeNeighborhood, myId, numProcs);
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read vector field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );

	// Create a vector field class from the file
	highF[0] = blockDim[0]-1.0f;
	highF[1] = blockDim[1]-1.0f;
	highF[2] = blockDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim, lowF, highF );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<VECTOR3>( vectorField, histogram );

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins at each point of the vector field ...\n" ); 
	starttime = ITL_util<float>::startTimer();	
	globalEntropyComputer->computeHistogramBinField( "vector", nBin );

	// Compute global entropy of the block
	globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	
	// Print block entropy (Process id = block id)
	printf( "Global entropy of block: %d, %f\n", myId, globalEntropyComputer->getGlobalEntropy() );

	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation Time: %f, %f seconds\n", myId, execTime[0], execTime[1] );
	else 			printf( "%d, %f, %f\n", myId, execTime[0], execTime[1] );
	
}// end function

void compute_localentropy_serial()
{
	// Read chunk of data from file
	if( verboseMode == 1 ) printf( "Reading vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinarySerial( vectorFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create a vector field class from the file
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim,
										   lowF, highF,
										   lowPad, highPad,
										   sizeNeighborhoodArray );

	// Initialize class that can compute entropy
	localEntropyComputer = new ITL_localentropy<VECTOR3>( vectorField, histogram );

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Entropy Computation
	if( verboseMode == 1 ) printf( "Computing entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeLocalEntropyOfField( nBin, true );
	execTime[2] = ITL_util<float>::endTimer( starttime );

	// Saving data to binary file
	if( verboseMode == 1 ) printf( "saving entropy field to binary file ...\n" );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	ITL_ioutil<float>::writeFieldBinarySerial( entropyField->getDataFull(), outFieldFile, entropyField->grid->dim, nDim );
	execTime[3] = ITL_util<float>::endTimer( starttime );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );
}// end function

void compute_localentropy_parallel()
{
	// Read relevant portion of data from file
	if( verboseMode == 1 ) printf( "%d: Reading vector field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel3( vectorFieldFile, nDim, dataDim,
												blockDim, nBlock,blockId,
												low, high,
												lowPad, highPad,
												sizeNeighborhoodArray, myId, numProcs );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create a vector field block from the relevant portion of the file
	for( int i=0; i<nDim; i++ )
	{
		lowF[i] = (float)low[i];
		highF[i] = (float)high[i];
	}
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim,
												  lowF, highF,
												  lowPad, highPad,
												  sizeNeighborhoodArray  );

	// Initialize class that can compute entropy
	localEntropyComputer = new ITL_localentropy<VECTOR3>( vectorField, histogram );

	if(verboseMode == 1 ) printf( "%d: Computing histogram at each point of the vector field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	
	if(verboseMode == 1 ) printf( "%d: Computing entropy at each point of the vector field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeLocalEntropyOfField( nBin, true );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	
	if(verboseMode == 1 ) printf( "%d: saving entropy field to binary file ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	ITL_ioutil<float>::writeFieldBinaryParallel2( entropyField->getDataFull(), outFieldFile, dataDim,
												blockDim, nBlock, blockId,
												lowF, nDim,
											        myId, numProcs );
	execTime[3] = ITL_util<float>::endTimer( starttime );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	if( verboseMode == 1 )  printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],execTime[2], execTime[3], execTime[4]  );
	else 			printf( "%d, %f, %f, %f\n", myId, execTime[0], execTime[4], execTime[3]  ); 	

}//end function

void compute_jointlocalentropy_serial()
{
	// Read chunk of data from both files 
	if( verboseMode == 1 ) printf( "Reading two vector fields ...\n" );
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinarySerial( vectorFieldFile, nDim, dataDim );
	data2 = ITL_ioutil<VECTOR3>::readFieldBinarySerial( vectorFieldFile2, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create first vector field class from the data
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim,
											lowF, highF,
											lowPad, highPad,
											sizeNeighborhood );
	vectorField2 = new ITL_field_regular<VECTOR3>( data2, nDim,
											lowF, highF,
											lowPad, highPad,
											sizeNeighborhood );

	// Initialize class that can compute entropy
	jointEntropyComputer = new ITL_localjointentropy<VECTOR3>( vectorField, vectorField2, histogram );

	// Histogram computation
	printf( "Converting vectors into joint histogram bins at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeJointHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Joint entropy Computation
	printf( "Computing joint entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeLocalJointEntropyOfField( nBin );
	execTime[2] = ITL_util<float>::endTimer( starttime );

	// Saving data to binary file
	printf( "saving entropy field to binary file ...\n" );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* jointEntropyField = jointEntropyComputer->getLocalJointEntropyField();
	ITL_ioutil<float>::writeFieldBinarySerial( jointEntropyField->getDataFull(), outFieldFile, jointEntropyField->grid->dim, nDim );
	execTime[3] = ITL_util<float>::endTimer( starttime );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	if( verboseMode == 1 )  printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],execTime[2], execTime[3], execTime[4]  );
	else 			printf( "%d, %f, %f, %f\n", myId, execTime[0], execTime[4], execTime[3]  ); 	

}// end function

void compute_jointlocalentropy_parallel()
{
	
	// Read chunk of data from first file
	if( verboseMode == 1 ) printf( "%d: Reading vector field ...\n", myId );
	if( myId == 0 )	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel3( vectorFieldFile, nDim, dataDim,
  														  blockDim, nBlock,blockId,
														  low, high,
														  lowPad, highPad,
														  sizeNeighborhoodArray, myId, numProcs );
	data2 = ITL_ioutil<VECTOR3>::readFieldBinaryParallel3( vectorFieldFile2, nDim, dataDim,
  														  blockDim, nBlock,blockId,
														  low, high,
														  lowPad, highPad,
														  sizeNeighborhoodArray, myId, numProcs );

	if( myId == 0 ) execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create first vector field class from the data
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim,
											lowF, highF,
											lowPad, highPad,
											sizeNeighborhoodArray );
	vectorField2 = new ITL_field_regular<VECTOR3>( data2, nDim,
											lowF, highF,
											lowPad, highPad,
											sizeNeighborhoodArray );

	// Initialize class that can compute entropy
	jointEntropyComputer = new ITL_localjointentropy<VECTOR3>( vectorField, vectorField2, histogram );

	// Histogram computation
	if( verboseMode == 1 ) printf( "%d: Converting vectors into joint histogram bins  ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeJointHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );

	// Joint entropy Computation
	if( verboseMode == 1 ) 	printf( "%d: Computing joint entropy at each point of the vector field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeLocalJointEntropyOfField( nBin );
	execTime[2] = ITL_util<float>::endTimer( starttime );

	// Saving data to binary file
	if( verboseMode == 1 ) printf( "%d: Saving joint local entropy field ...\n", myId );
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	ITL_ioutil<float>::writeFieldBinaryParallel2( jointEntropyComputer->getLocalJointEntropyField()->getDataFull(),
											  	  (char*)outFieldFile, dataDim,
												  blockDim, nBlock, blockId,
												  lowF, nDim,
												  myId, numProcs );
	if( myId == 0 ) execTime[3] = ITL_util<float>::endTimer( starttime );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );
}//end function

