/**
 * @file MainIT_regscalar.cpp
 * Application program for entropy computation of regular scalar field.
 * Created on: April 28, 2011
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
#include "ITL_globaljointentropy.h"
#include "ITL_field_regular.h"

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

int numProcs;
int myId;

int functionType = 0;
int nDim = 3;
int nBin = 1000;
int dataDim[3];
int low[3];
int high[3];
float lowF[] = {0, 0, 0};
float highF[] = {0, 0, 0};
int lowPad[] = {0, 0, 0};
int highPad[] = {0, 0, 0};
int sizeNeighborhood = 0;
int sizeNeighborhoodArray[] = {0, 0, 0};
int blockDim[3];
int nBlock[3];
int blockId[3];
int method = 0;
int verboseMode = 1;
float histogramLowEnd = 0;
float histogramHighEnd = 0;

SCALAR *scalarFieldData = NULL;
SCALAR *scalarFieldData2 = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<SCALAR> *scalarField2 = NULL;

ITL_localentropy<SCALAR> *localEntropyComputer = NULL;
ITL_globalentropy<SCALAR> *globalEntropyComputer = NULL;
ITL_globaljointentropy<SCALAR> *globalJointEntropyComputer = NULL;

const char *scalarFieldFile = NULL;
const char *scalarFieldFile2 = NULL;
const char *heightFieldFile = NULL;
const char *tetMeshFile = NULL;
const char *outFieldFile = NULL;

double execTime[5];
clock_t starttime, endtime;

/**
 * Serial global entropy computation function for regular scalar field.
 */
void compute_globalentropy_serial();
/**
 * Parallel global entropy computation function for regular scalar field.
 */
void compute_globalentropy_parallel();
/**
 * Parallel blockwise global entropy computation function for regular scalar field.
 */
void compute_blockwiseglobalentropy_parallel();
/**
 * Serial local entropy computation function for regular scalar field.
 */
void compute_localentropy_regularfield_serial();
/**
 * Parallel local entropy computation function for regular scalar field.
 */
void compute_localentropy_parallel();
/**
 * Serial global joint entropy computation function for regular scalar field.
 */
void compute_globaljointentropy_regularfield_serial();

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
	scalarFieldFile =  ITL_util<float>::getArgWithName( "scalarField", &argNames, &argValues );
	scalarFieldFile2 =  ITL_util<float>::getArgWithName( "scalarField2", &argNames, &argValues );
	outFieldFile =  ITL_util<float>::getArgWithName( "outField", &argNames, &argValues );
	nDim = atoi(  ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	nBin = atoi(  ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
	sizeNeighborhood = atoi( ITL_util<float>::getArgWithName( "neighborhoodSize", &argNames, &argValues ) );
	sizeNeighborhoodArray[0] = atoi( ITL_util<float>::getArgWithName( "neighborhoodSizeX", &argNames, &argValues ) );
	sizeNeighborhoodArray[1] = atoi(  ITL_util<float>::getArgWithName( "neighborhoodSizeY", &argNames, &argValues ) );
	sizeNeighborhoodArray[2] = atoi(  ITL_util<float>::getArgWithName( "neighborhoodSizeZ", &argNames, &argValues ) );
	nBlock[0] = atoi(  ITL_util<float>::getArgWithName( "nBlockX", &argNames, &argValues ) );
	nBlock[1] = atoi(  ITL_util<float>::getArgWithName( "nBlockY", &argNames, &argValues ) );
	nBlock[2] = atoi(  ITL_util<float>::getArgWithName( "nBlockZ", &argNames, &argValues ) );
	method = atoi(  ITL_util<float>::getArgWithName( "method", &argNames, &argValues ) );
	double d =  atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramLowEnd = (float)atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramHighEnd = (float)atof( ITL_util<float>::getArgWithName( "histHigh", &argNames, &argValues ) );
	verboseMode = atoi( ITL_util<float>::getArgWithName( "verbose", &argNames, &argValues ) );

	switch( functionType )
	{
	case 0:
		printf( "Entering serial computation of global entropy of the regular scalar field ...\n" );	
		compute_globalentropy_serial();
		break;
	case 1:
		printf( "Entering parallel computation of global entropy of the regular scalar field ...\n" );	
		compute_globalentropy_parallel();
		break;
	case 2:
		printf( "Entering parallel computation of blockwise global entropy field from the regular scalar field ...\n" );	
		compute_blockwiseglobalentropy_parallel();
		break;
	case 3:
		printf( "Entering serial computation of local entropy field from the regular scalar field ...\n" );	
		compute_localentropy_regularfield_serial();
		break;
	case 4:
		printf( "Entering parallel computation of local entropy field from the regular scalar field ...\n" );	
		compute_localentropy_parallel();
		break;	
	case 5:
		printf( "Entering serial computation of global joint entropy of the regular scalar field ...\n" );	
		compute_globaljointentropy_regularfield_serial();
		break;
	default:
		break;
	}// end switch

	// Clear up
	if( localEntropyComputer != NULL ) delete localEntropyComputer;
	if( globalEntropyComputer != NULL ) delete globalEntropyComputer;
	//if( scalarField != NULL ) delete scalarField;
		
	// Finalize MPI
	MPI_Finalize();

}// end main

void compute_globalentropy_serial()
{	
	// Read scalar field
	if( verboseMode == 1 ) printf( "Reading scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<float>::readFieldBinarySerial( scalarFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read scalar field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );
	
	// Create a scalar field class from the file
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF, lowPad, highPad, sizeNeighborhoodArray );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<SCALAR>( scalarField );

	// Global entropy computation using different methods
	// 0: Histogram based
	// 1: KDE based
	if( method == 0 )
	{
		if( histogramLowEnd != histogramHighEnd )
			globalEntropyComputer->setHistogramRange( histogramLowEnd, histogramHighEnd );	
		
		// Cross-validation
		//globalEntropyComputer->crossValidate( "scalar", 1000, 100, 10 );
		
		
		// Histogram computation
		if( verboseMode == 1 ) printf( "Converting scalars into histogram bins ...\n" );
		starttime = ITL_util<float>::startTimer();
		globalEntropyComputer->computeHistogramBinField( "scalar", nBin );
		execTime[1] = ITL_util<float>::endTimer( starttime );
		
		// Entropy Computation
		if( verboseMode == 1 ) printf( "Computing global entropy of the scalar field ...\n" );
		globalEntropyComputer->computeGlobalEntropyOfField( nBin, false );

		// Get histogram frequencies
		int *freqList = new int[nBin];
		globalEntropyComputer->getHistogramFrequencies( nBin, freqList );
		if( verboseMode == 1 ) 
		{
			for( int p=0; p<nBin; p++ )
				printf( "f[%d] = %d, ", p, freqList[p] );
			printf( "\n" );
		}

		float globalEntropy = globalEntropyComputer->getGlobalEntropy();
		execTime[1] = ITL_util<float>::endTimer( starttime );
		printf( "Global Entropy: %f\n", globalEntropy );
	}
	else if( method == 1 )
	{
		// KDE based entropy computation
		printf( "Computing global entropy of the scalar field based on KDE...\n" );
		starttime = ITL_util<float>::startTimer();
		globalEntropyComputer->computeGlobalEntropyOfField( nBin, true, 1 );
		float globalEntropy = globalEntropyComputer->getGlobalEntropy();
		execTime[1] = ITL_util<float>::endTimer( starttime );
		printf( "Global Entropy: %f\n", globalEntropy );
	}

	// Runtime
	if( verboseMode == 1 )  printf( "%d: Read/Total Computation Time: %f, %f seconds\n", myId, execTime[0], execTime[1] );
	else			printf( "%d, %f, %f\n", myId, execTime[0], execTime[1] );	

	delete scalarField;
	
}// end function

void compute_globalentropy_parallel()
{
	// Read chunk of data from file
	if( verboseMode == 1 ) printf( "Reading scalar field ...\n" );	
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinaryParallel2( scalarFieldFile, nDim, dataDim,
								blockDim, nBlock,blockId,
								low, high,
								lowPad, highPad,
								sizeNeighborhood, myId, numProcs);
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read scalar field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );

	// Create a vector field class from the file
	highF[0] = blockDim[0]-1.0f;
	highF[1] = blockDim[1]-1.0f;
	highF[2] = blockDim[2]-1.0f;
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<SCALAR>( scalarField );

	if( histogramLowEnd != histogramHighEnd )
		globalEntropyComputer->setHistogramRange( histogramLowEnd, histogramHighEnd );				

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins at each point of the vector field ...\n" ); 
	starttime = ITL_util<float>::startTimer();	
	globalEntropyComputer->computeHistogramBinField( "scalar", nBin );
	if( verboseMode == 1 ) printf( "Done\n" );

	// Compute frequencies
	globalEntropyComputer->computeHistogramFrequencies( nBin );
	
	// Get histogram frequencies
	#if defined( _WIN32 ) || defined( _WIN64 )
		int* freqList = new int[nBin];
	#endif
	#if defined( _LINUX )
		int freqList[nBin];
	#endif

	globalEntropyComputer->getHistogramFrequencies( nBin, freqList );
	if( verboseMode == 1 ) 
	{	
		for( int i=0; i<nBin; i++ )
			printf( "%d:f[%d]=%d\n", myId, i, freqList[i] );
	}

	// Sync with all processors and and sum up all histogram frequencies to processor 0
	MPI_Barrier( MPI_COMM_WORLD );

	#if defined( _WIN32 ) || defined( _WIN64 )
		int* reducedFreqList = new int[nBin];
	#endif
	#if defined( _LINUX )
		int reducedFreqList[nBin];
	#endif
	MPI_Reduce( freqList, reducedFreqList, nBin, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );  

	// Compute and print global entropy at processor 0	
	if( myId == 0 )	
	{
		float globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( reducedFreqList, dataDim[0]*dataDim[1]*dataDim[2], nBin, false );
		printf( "Global entropy of scalar field: %f\n", globalEntropy );
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
	if( verboseMode == 1 ) printf( "Reading scalar field ...\n" );	
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinaryParallel2( scalarFieldFile, nDim, dataDim,
								blockDim, nBlock,blockId,
								low, high,
								lowPad, highPad,
								sizeNeighborhood, myId, numProcs);
	execTime[0] = ITL_util<float>::endTimer( starttime );
	if( verboseMode == 1 ) printf( "Read scalar field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );

	// Create a vector field class from the file
	highF[0] = blockDim[0]-1.0f;
	highF[1] = blockDim[1]-1.0f;
	highF[2] = blockDim[2]-1.0f;
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<SCALAR>( scalarField );

	if( histogramLowEnd != histogramHighEnd )
		globalEntropyComputer->setHistogramRange( histogramLowEnd, histogramHighEnd );				

	// Histogram computation
	if( verboseMode == 1 ) printf( "Converting vectors into histogram bins at each point of the vector field ...\n" ); 
	starttime = ITL_util<float>::startTimer();	
	globalEntropyComputer->computeHistogramBinField( "scalar", nBin );
	if( verboseMode == 1 ) printf( "Done\n" );

	// Compute global entropy of the block
	globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	
	// Print block entropy (Process id = block id)
	printf( "Global entropy of block: %d, %f\n", myId, globalEntropyComputer->getGlobalEntropy() );

	// Runtime
	if( verboseMode == 1 ) 	printf( "%d: Read/Computation Time: %f, %f seconds\n", myId, execTime[0], execTime[1] );
	else 			printf( "%d, %f, %f\n", myId, execTime[0], execTime[1] );
	
}// end function


void compute_localentropy_regularfield_serial()
{
	// Initialize ITL
	ITL_base::ITL_init();

	// Create a scalar field class from the file
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;

	// Read scalar field
	printf( "Reading scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<float>::readFieldBinarySerial( scalarFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );
	printf( "Read scalar field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );


	// Initialize scalar field class
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim,
 							lowF, highF,
							lowPad, highPad,
							sizeNeighborhoodArray );

	// Initialize class that can compute local entropy field
	localEntropyComputer = new ITL_localentropy<SCALAR>( scalarField );

	// Histogram computation
	printf( "Converting scalars into histogram bins at each point of the scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "scalar", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Entropy Computation
	printf( "Computing entropy at each point of the scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeEntropyOfField( nBin, false );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Saving data to binary file
	printf( "Saving local entropy field to binary file ...\n" );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	ITL_ioutil<float>::writeFieldBinarySerial( entropyField->getDataFull(), outFieldFile, entropyField->grid->dim, nDim );
	execTime[3] = ITL_util<float>::endTimer( starttime );
	//printf( "Done\n" );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );

	delete scalarField;
	delete entropyField;	

}// end function

void compute_localentropy_parallel()
{
	// Read relevant portion of data from file
	if(verboseMode == 1 ) printf( "%d:Reading portion of scalar field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	//scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinaryParallel2( scalarFieldFile, nDim, dataDim, blockDim, nBlock,blockId, low, high, lowPad, highPad, sizeNeighborhood, myId, numProcs );
	scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinaryParallel3( scalarFieldFile, nDim, dataDim, blockDim, nBlock,blockId, low, high, lowPad, highPad, sizeNeighborhoodArray, myId, numProcs );
	execTime[0] = ITL_util<float>::endTimer( starttime );
	for( int i=0; i<blockDim[0]*blockDim[1]*blockDim[2]; i++ )
	{
		if( scalarFieldData[i] == HUGE_VAL )
			scalarFieldData[i] = 0.0f;
	}

	// Create a scalar field block from the relevant portion of the file
	for( int i=0; i<nDim; i++ )
	{
		lowF[i] = (float)low[i];
		highF[i] = (float)high[i];
	}
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF, lowPad, highPad, sizeNeighborhoodArray  );

	// Initialize class that can compute entropy
	localEntropyComputer = new ITL_localentropy<SCALAR>( scalarField );

	if(verboseMode == 1 )	printf( "%d: Computing histogram at each point of the scalar field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "scalar", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );

	if(verboseMode == 1 ) printf( "%d: Computing entropy at each point of the scalar field ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeEntropyOfField( nBin, true );
	execTime[2] = ITL_util<float>::endTimer( starttime );


	if(verboseMode == 1 ) printf( "%d: saving entropy field to binary file ...\n", myId );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	ITL_ioutil<float>::writeFieldBinaryParallel2( entropyField->getDataFull(), outFieldFile, dataDim, blockDim, nBlock, blockId, lowF, nDim, myId, numProcs );
	execTime[3] = ITL_util<float>::endTimer( starttime );

	execTime[4] = execTime[1] + execTime[2];
	if(verboseMode == 1 ) printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1], execTime[2], execTime[3], execTime[4] );
	else		      printf( "%d, %f, %f, %f\n", myId, execTime[0], execTime[4], execTime[3] );						 				

	delete entropyField;
	delete scalarField;

}//end function



// ADD-BY-Abon 07/19/2011-BEGIN
void compute_globaljointentropy_regularfield_serial()
{
	// Read chunk of data from both files 
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinarySerial( scalarFieldFile, nDim, dataDim );
	scalarFieldData2 = ITL_ioutil<SCALAR>::readFieldBinarySerial( scalarFieldFile2, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create two scalar field classes from the data
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF, lowPad, highPad, sizeNeighborhood );
	scalarField2 = new ITL_field_regular<SCALAR>( scalarFieldData2, nDim, lowF, highF, lowPad, highPad, sizeNeighborhood );

	// Initialize class that can compute global joint entropy
	globalJointEntropyComputer = new ITL_globaljointentropy<SCALAR>( scalarField, scalarField2 );

	// Histogram computation
	printf( "Converting scalars into joint histogram bins at each point of the scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	if( histogramLowEnd != histogramHighEnd )
		globalJointEntropyComputer->setJointHistogramRange( histogramLowEnd, histogramHighEnd, histogramLowEnd, histogramHighEnd );	
	globalJointEntropyComputer->computeJointHistogramBinField( "scalar", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Global Joint entropy Computation
	printf( "Computing joint entropy of the entire scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	globalJointEntropyComputer->computeGlobalJointEntropyOfField( nBin, false );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );
	
	// Display computed global joint entropy
	printf( "Joint global entropy of scalar fields: %f\n", globalJointEntropyComputer->getGlobalJointEntropy() );
	

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1], execTime[2], execTime[3], execTime[4]  );

}// end function
// ADD-BY-Abon 07/19/2011-END




