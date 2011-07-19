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

ITL_field_regular<VECTOR3> *vectorField = NULL;
ITL_field_regular<VECTOR3> *vectorField2 = NULL;

ITL_localentropy<VECTOR3> *localEntropyComputer = NULL;
ITL_globalentropy<VECTOR3> *globalEntropyComputer = NULL;
ITL_localjointentropy<VECTOR3> *jointEntropyComputer = NULL;

/**
 * Serial converter from vec to raw format
 * @param nArg narg
 * @param argV argv
 */
void create_raw_file( int nArg, char** argV );
/**
 * Serial entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_localentropy_sequential();
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
void compute_jointlocalentropy_sequential();
/**
 * Parallel joint entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_jointlocalentropy_parallel();
/**
 * Serial global entropy computation function with runtime computation.
 * @param nArg narg
 * @param argV argv
 */
void compute_globalentropy_sequential();

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
	
	// Initialize ITL
	ITL_base::ITL_init();

	// Initialize histogram
	ITL_histogram::ITL_init_histogram( patchFile, nBin );
		
	switch( functionType )
	{
	case 0:
		create_raw_file( argc, argv );
		break;
	case 1:
		printf( "Entering serial computation of local entropy field from the regular vector field ...\n" );	
		compute_localentropy_sequential();
		break;
	case 2:
		printf( "Entering parallel computation of local entropy field from the regular vector field ...\n" );	
		compute_localentropy_parallel();
		break;
	case 3:
		printf( "Entering serial computation of local joint entropy field from two regular vector fields ...\n" );	
		compute_jointlocalentropy_sequential();
		break;
	case 4:
		printf( "Entering parallel computation of local joint entropy fields from two regular vector fields ...\n" );	
		compute_jointlocalentropy_parallel();
		break;
	case 5:
		printf( "Entering serial computation of global entropy field from a regular vector field ...\n" );
		compute_globalentropy_sequential();
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

void compute_localentropy_sequential()
{
	// Read chunk of data from file
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
	localEntropyComputer = new ITL_localentropy<VECTOR3>( vectorField );

	// Histogram computation
	printf( "Converting vectors into histogram bins at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Entropy Computation
	printf( "Computing entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeEntropyOfField( nBin, true );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Saving data to binary file
	printf( "saving entropy field to binary file ...\n" );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* entropyField = localEntropyComputer->getEntropyField();
	cout << entropyField->grid->dim[0] << " " << entropyField->grid->dim[1] << " " << entropyField->grid->dim[2] << endl;
	
	ITL_ioutil<float>::writeFieldBinarySerial( entropyField->getDataFull(), outFieldFile, entropyField->grid->dim, nDim );
	execTime[3] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );
}// end function

void compute_globalentropy_sequential()
{
	// Read chunk of data from file
	starttime = ITL_util<float>::startTimer();
	data = ITL_ioutil<VECTOR3>::readFieldBinarySerial( vectorFieldFile, nDim, dataDim );
	execTime[0] = ITL_util<float>::endTimer( starttime );

	// Create a vector field class from the file
	highF[0] = dataDim[0]-1.0f;
	highF[1] = dataDim[1]-1.0f;
	highF[2] = dataDim[2]-1.0f;
	vectorField = new ITL_field_regular<VECTOR3>( data, nDim, lowF, highF, lowPad, highPad, sizeNeighborhoodArray );

	// Initialize class that can compute entropy
	globalEntropyComputer = new ITL_globalentropy<VECTOR3>( vectorField );

	// Histogram computation
	printf( "Converting vectors into histogram bins at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	globalEntropyComputer->computeHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Global entropy Computation
	printf( "Computing entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	globalEntropyComputer->computeGlobalEntropyOfField( nBin, false );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Print global entropy
	printf( "Global entropy of vector field: %f\n", globalEntropyComputer->getGlobalEntropy() );

	// Runtime
	execTime[3] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Total Computation Time: %f %f %f %f seconds\n", myId, execTime[0], execTime[1], execTime[2], execTime[3] );
	
}// end function

void compute_localentropy_parallel()
{
	// Read relevant portion of data from file
	if( myId == 0 )	starttime = ITL_util<float>::startTimer();
	//data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel2( vectorFieldFile, nDim, dataDim,
  	//												blockDim, nBlock,blockId,
	//												low, high,
	//												lowPad, highPad,
	//												sizeNeighborhood, myId, numProcs );
	data = ITL_ioutil<VECTOR3>::readFieldBinaryParallel3( vectorFieldFile, nDim, dataDim,
												blockDim, nBlock,blockId,
												low, high,
												lowPad, highPad,
												sizeNeighborhoodArray, myId, numProcs );

	if( myId == 0 ) execTime[0] = ITL_util<float>::endTimer( starttime );

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
	localEntropyComputer = new ITL_localentropy<VECTOR3>( vectorField );

	printf( "%d: Computing histogram at each point of the vector field ...\n", myId );
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	localEntropyComputer->computeHistogramBinField( "vector", nBin );
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
	//ITL_ioutil<float>::writeFieldBinaryParallel( entropyField->getDataFull(), outFieldFile,
	//										     dataDim, blockDim,
	//										     lowF, nDim,
	//										     myId, numProcs );
	ITL_ioutil<float>::writeFieldBinaryParallel2( entropyField->getDataFull(), outFieldFile, dataDim,
												blockDim, nBlock, blockId,
												lowF, nDim,
											        myId, numProcs );
	if( myId == 0 ) execTime[3] = ITL_util<float>::endTimer( starttime );
	//printf( "%d: Done\n", myId );

	execTime[4] = execTime[1] + execTime[2];
	if( myId == 0 ) printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n",
																myId, execTime[0], execTime[1],
																execTime[2], execTime[3], execTime[4] );
}//end function

void compute_jointlocalentropy_sequential()
{
	// Read chunk of data from both files 
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
	jointEntropyComputer = new ITL_localjointentropy<VECTOR3>( vectorField, vectorField2 );

	// Histogram computation
	printf( "Converting vectors into joint histogram bins at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeJointHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Joint entropy Computation
	printf( "Computing joint entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeLocalJointEntropyOfField( nBin );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Saving data to binary file
	printf( "saving entropy field to binary file ...\n" );
	starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* jointEntropyField = jointEntropyComputer->getLocalJointEntropyField();
	ITL_ioutil<float>::writeFieldBinarySerial( jointEntropyField->getDataFull(), outFieldFile, jointEntropyField->grid->dim, nDim );
	execTime[3] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );

}// end function

void compute_jointlocalentropy_parallel()
{
	
	// Read chunk of data from first file
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
	jointEntropyComputer = new ITL_localjointentropy<VECTOR3>( vectorField, vectorField2 );

	// Histogram computation
	printf( "Converting vectors into joint histogram bins at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeJointHistogramBinField( "vector", nBin );
	execTime[1] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Joint entropy Computation
	printf( "Computing joint entropy at each point of the vector field ...\n" );
	starttime = ITL_util<float>::startTimer();
	jointEntropyComputer->computeLocalJointEntropyOfField( nBin );
	execTime[2] = ITL_util<float>::endTimer( starttime );
	printf( "Done\n" );

	// Saving data to binary file
	if( myId == 0 ) starttime = ITL_util<float>::startTimer();
	ITL_field_regular<float>* jointEntropyField = jointEntropyComputer->getLocalJointEntropyField();
	//ITL_ioutil<float>::writeFieldBinaryParallel( entropyField->getDataFull(), outFieldFile,
	//										     dataDim, blockDim,
	//										     lowF, nDim,
	//										     myId, numProcs );

	ITL_ioutil<float>::writeFieldBinaryParallel2( jointEntropyField->getDataFull(), (char*)outFieldFile, dataDim,
													blockDim, nBlock, blockId,
													lowF, nDim,
													myId, numProcs );													
	if( myId == 0 ) execTime[3] = ITL_util<float>::endTimer( starttime );
	printf( "%d: Done\n", myId );

	// Runtime
	execTime[4] = execTime[1] + execTime[2];
	printf( "%d: Read/Histogram/Entropy/Write/Total Computation Time: %f %f %f %f %f seconds\n", myId, execTime[0], execTime[1],
																			execTime[2], execTime[3], execTime[4]  );
}//end function

