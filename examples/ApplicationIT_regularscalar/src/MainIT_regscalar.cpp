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

SCALAR *scalarFieldData = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_localentropy<SCALAR> *localEntropyComputer = NULL;
ITL_globalentropy<SCALAR> *globalEntropyComputer = NULL;

const char *scalarFieldFile = NULL;
const char *heightFieldFile = NULL;
const char *tetMeshFile = NULL;
const char *outFieldFile = NULL;

double execTime[5];
clock_t starttime, endtime;

void getArgs( const char *argsFileName );
const char* getArgWithName( const char* name );

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
	getArgs( argv[1] );

	// Parse command line arguments
	functionType = atoi( getArgWithName( "functionType" ) );
	scalarFieldFile = getArgWithName( "scalarField" );
	outFieldFile = getArgWithName( "outField" );
	nDim = atoi( getArgWithName( "nDim" ) );
	nBin = atoi( getArgWithName( "nBin" ) );
	sizeNeighborhood = atoi( getArgWithName( "neighborhoodSize" ) );
	sizeNeighborhoodArray[0] = atoi( getArgWithName( "neighborhoodSizeX" ) );
	sizeNeighborhoodArray[1] = atoi( getArgWithName( "neighborhoodSizeY" ) );
	sizeNeighborhoodArray[2] = atoi( getArgWithName( "neighborhoodSizeZ" ) );
	nBlock[0] = atoi( getArgWithName( "nBlockX" ) );
	nBlock[1] = atoi( getArgWithName( "nBlockY" ) );
	nBlock[2] = atoi( getArgWithName( "nBlockZ" ) );
	method = atoi( getArgWithName( "method" ) );

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
	default:
		break;
	}// end switch

	// Clear up
	if( scalarFieldData != NULL ) delete [] scalarFieldData;
		
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

void getArgs( const char *argsFileName )
{
	char argName[100];
	char argVal[200];

	// Open file containing list of arguments
	FILE* argsFile = fopen( argsFileName, "r" );

	// Scan the file to create two lists of strings
	// List of argument names and list of argument values
	while( true )
	{
		// Read two strings in the current line
		fscanf( argsFile, "%s %s", argName, argVal );

		// Break if end of file reached
		if( strcmp( argName, "EOF") == 0 )
			break;

		// Place strings into corresponding lists
		string name( argName );
		string val( argVal );

		// Push strings in to lists
		argNames.push_back( name );
		argValues.push_back( val );

	}

	// Close file
	fclose( argsFile );

}// end function

const char* getArgWithName( const char* name )
{
	list<string>::iterator iterName;
	list<string>::iterator iterVal;

	for( iterName = argNames.begin(), iterVal = argValues.begin();
		 iterName != argNames.end();
		 ++iterName, ++iterVal )
	{
		if( strcmp( (*iterName).c_str(), name ) == 0 )
			return (*iterVal).c_str();

	}// end for

	return "";

}// end function



