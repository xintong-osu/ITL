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
int sizeNeighborhood = 6;
int sizeNeighborhoodArray[] = {0, 0, 0};
int blockDim[3];
int nBlock[3];
int blockId[3];
int method = 0;
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
 * Serial local entropy computation function for regular scalar field.
 */
void compute_localentropy_regularfield_serial();
/**
 * Parallel local entropy computation function for regular scalar field.
 */
void compute_localentropy_regularfield_parallel();
/**
 * Serial global entropy computation function for regular scalar field.
 */
void compute_globalentropy_regularfield_serial();
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
	histogramLowEnd = atof( ITL_util<float>::getArgWithName( "histLow", &argNames, &argValues ) );
	histogramHighEnd = atof( ITL_util<float>::getArgWithName( "histHigh", &argNames, &argValues ) );

	switch( functionType )
	{
	case 0:
		printf( "Entering serial computation of local entropy field from the regular scalar field ...\n" );	
		compute_localentropy_regularfield_serial();
		break;
	case 1:
		printf( "Entering parallel computation of local entropy field from the regular scalar field ...\n" );	
		compute_localentropy_regularfield_parallel();
		break;
	case 2:
		printf( "Entering serial computation of global entropy of the regular scalar field ...\n" );	
		compute_globalentropy_regularfield_serial();
		break;
	// ADD-BY-Abon 07/19/2011-BEGIN	
	case 3:
		printf( "Entering serial computation of global joint entropy of the regular scalar field ...\n" );	
		compute_globaljointentropy_regularfield_serial();
		break;
	// ADD-BY-Abon 07/19/2011-END
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
	localEntropyComputer->computeEntropyOfField( nBin, true );
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

void compute_localentropy_regularfield_parallel()
{
	// Read relevant portion of data from file
	printf( "%d:Reading portion of scalar field ...\n", myId );
	if( myId == 0 )	starttime = ITL_util<float>::startTimer();
	//data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel2( vectorFieldFile, nDim, dataDim, blockDim, nBlock,blockId, low, high, lowPad, highPad, sizeNeighborhood, myId, numProcs );
	scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinaryParallel3( scalarFieldFile, nDim, dataDim, blockDim, nBlock,blockId, low, high, lowPad, highPad, sizeNeighborhoodArray, myId, numProcs );
	if( myId == 0 ) execTime[0] = ITL_util<float>::endTimer( starttime );
	printf( "Done.\n" );

	for( int i=0; i<blockDim[0]*blockDim[1]*blockDim[2]; i++ )
	{
		if( scalarFieldData[i] == HUGE_VAL )
			scalarFieldData[i] = 0.0f;
	}
	printf("replaced\n");

	// Create a scalar field block from the relevant portion of the file
	for( int i=0; i<nDim; i++ )
	{
		lowF[i] = (float)low[i];
		highF[i] = (float)high[i];
	}
	scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF, lowPad, highPad, sizeNeighborhoodArray  );

	// Initialize class that can compute entropy
	localEntropyComputer = new ITL_localentropy<SCALAR>( scalarField );

	printf( "%d: Computing histogram at each point of the vector field ...\n", myId );
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "scalar", nBin );
	if( myId == 0 ) execTime[1] = ITL_util<float>::endTimer( starttime );
	//printf( "%d: Done\n", myId );

	printf( "%d: Computing entropy at each point of the vector field ...\n", myId );
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeEntropyOfField( nBin, true );
	if( myId == 0 ) execTime[2] = ITL_util<float>::endTimer( starttime );
	//printf( "%d: Done\n", myId );

	printf( "%d: saving entropy field to binary file ...\n", myId );
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	//ITL_ioutil<float>::writeFieldBinaryParallel( entropyField->getDataFull(), outFieldFile, dataDim, blockDim, lowF, nDim, myId, numProcs );
	ITL_ioutil<float>::writeFieldBinaryParallel2( entropyField->getDataFull(), outFieldFile, dataDim, blockDim, nBlock, blockId, lowF, nDim, myId, numProcs );
	if( myId == 0 ) execTime[3] = ITL_util<float>::endTimer( starttime );
	//printf( "%d: Done\n", myId );

	execTime[4] = execTime[1] + execTime[2];
	if( myId == 0 ) printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1], execTime[2], execTime[3], execTime[4] );

	delete entropyField;
	delete scalarField;

}//end function

void compute_globalentropy_regularfield_serial()
{	
	// Read scalar field
	printf( "Reading scalar field ...\n" );
	starttime = ITL_util<float>::startTimer();
	scalarFieldData = ITL_ioutil<float>::readFieldBinarySerial( scalarFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );
	printf( "Read scalar field of dimension: %d %d %d\n", dataDim[0], dataDim[1], dataDim[2] );
	
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
		// Histogram computation
		printf( "Converting scalars into histogram bins ...\n" );
		starttime = ITL_util<float>::startTimer();
		globalEntropyComputer->computeHistogramBinField( "scalar", nBin );
		execTime[1] = ITL_util<float>::endTimer( starttime );
		printf( "Done\n" );

		// Entropy Computation
		printf( "Computing global entropy of the scalar field ...\n" );
		starttime = ITL_util<float>::startTimer();
		globalEntropyComputer->computeGlobalEntropyOfField( nBin, false );
		float globalEntropy = globalEntropyComputer->getGlobalEntropy();
		execTime[2] = ITL_util<float>::endTimer( starttime );
		printf( "Done\nGlobal Entropy: %f\n", globalEntropy );
	}
	else if( method == 1 )
	{
		// KDE based entropy computation
		printf( "Computing global entropy of the scalar field based on KDE...\n" );
		starttime = ITL_util<float>::startTimer();
		globalEntropyComputer->computeGlobalEntropyOfField( nBin, true, 1 );
		float globalEntropy = globalEntropyComputer->getGlobalEntropy();
		execTime[2] = ITL_util<float>::endTimer( starttime );
		printf( "Done\nGlobal Entropy: %f\n", globalEntropy );
	}

	// Runtime
	execTime[3] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Total Computation Time: %f %f %f %f seconds\n", myId, execTime[0], execTime[1], execTime[2], execTime[3]  );

	delete scalarField;
	
}// end function

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




