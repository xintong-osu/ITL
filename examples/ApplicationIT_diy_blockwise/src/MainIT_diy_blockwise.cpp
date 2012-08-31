/**
 * @file MainIT_diy.cpp
 * Application program for entropy computation using diy framework
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
#include "ITL_globalentropy.h"

#define NO_DIY 0
#define WITH_DIY 1 

using namespace std;

// Global variables
list<string> argNames;
list<string> argValues;

char* inputFieldFile = NULL; 	
char* outFile = NULL;
char* patchFile = NULL;
SCALAR *globalEntropyList = NULL;

SCALAR *scalarFieldData = NULL;
VECTOR3 *vectorFieldData = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR>* histMapper_scalar = NULL;
ITL_histogrammapper<VECTOR3>* histMapper_vector = NULL;

ITL_field_regular<SCALAR> *scalarField = NULL;
ITL_field_regular<VECTOR3> *vectorField = NULL;

ITL_field_regular<int> *binField = NULL;

ITL_globalentropy<SCALAR> *globalEntropyComputer_scalar = NULL;
ITL_globalentropy<VECTOR3> *globalEntropyComputer_vector = NULL;

ITL_field_regular<SCALAR> *subScalarFieldArray = NULL;
ITL_field_regular<VECTOR3> *subVectorFieldArray = NULL;


int tot_blocks = 512;
int nblocks = 0;				// My local number of blocks

double execTime[5];
clock_t starttime, endtime;

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
void *CreateType(void *item, DIY_Datatype *dtype)
{
	//MPI_Datatype *dtype = new MPI_Datatype; // DIY will free this resource for you
	//MPI_Type_contiguous(nBin, MPI_INT, dtype);
	//MPI_Type_commit(dtype); // DIY will free this resource for you
	//abs_addr = false;
	//return dtype;
	//struct map_block_t map[1] = {
	//	{MPI_FLOAT, OFST, nblocks, 0, 1},
	//};
	//DIY_Create_datatype(DIY_Addr(item), 1, map, dtype);
	//return MPI_BOTTOM;
  
	DIY_Create_vector_datatype( nblocks, 1, DIY_FLOAT, dtype );
	return item;
}

/**
 * Main function.
 * Program starts from here.
 */
int main( int argc, char** argv )
{
	int numProcs;
	int rank;
	int num_threads = 4; // number of threads DIY can use

	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	int nDim;
	int nBin;
	
	int dataSize[3];

	int fieldType = -1;
	int method = 0;
	int parallelOpMode = NO_DIY;
	int verboseMode = 1;

	float lowF[3];
	float highF[3];
	int nPartition[3]; 
	float histogramLowEnd = 0;
	float histogramHighEnd = 0;

	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	inputFieldFile = ITL_util<float>::getArgWithName( "inputField", &argNames, &argValues );
	outFile =  ITL_util<float>::getArgWithName( "outFile", &argNames, &argValues );
	patchFile =  ITL_util<float>::getArgWithName( "patchFile", &argNames, &argValues );
	fieldType = atoi( ITL_util<float>::getArgWithName( "fieldType", &argNames, &argValues ) );
	nDim = atoi( ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	nBin = atoi( ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
	parallelOpMode = atoi( ITL_util<float>::getArgWithName( "useDIY", &argNames, &argValues ) );
	dataSize[0] = atoi( ITL_util<float>::getArgWithName( "nX", &argNames, &argValues ) );
	dataSize[1] = atoi( ITL_util<float>::getArgWithName( "nY", &argNames, &argValues ) );
	dataSize[2] = atoi( ITL_util<float>::getArgWithName( "nZ", &argNames, &argValues ) );
	tot_blocks = atoi( ITL_util<float>::getArgWithName( "nBlock", &argNames, &argValues ) );
	nPartition[0] = atoi( ITL_util<float>::getArgWithName( "nPartX", &argNames, &argValues ) );
	nPartition[1] = atoi( ITL_util<float>::getArgWithName( "nPartY", &argNames, &argValues ) );
	nPartition[2] = atoi( ITL_util<float>::getArgWithName( "nPartZ", &argNames, &argValues ) );
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
	//cout << "-5" << endl;

	// Open file to dump output only if program running on single processor
	FILE *dumpFile = NULL;	
	if( numProcs == 1 && verboseMode == 1 )	dumpFile = fopen( "./dump.csv", "w" );
	//cout << "-4" << endl;

	if( parallelOpMode == NO_DIY )
	{
		lowF[0] = lowF[1] = lowF[2] = 0;
		highF[0] = (float)(dataSize[0]-1);
		highF[1] = (float)(dataSize[1]-1);
		highF[2] = (float)(dataSize[2]-1);
		//cout << "-3" << endl;

		// Compute total number of blocks to be created
		nblocks = nPartition[0] * nPartition[1] * nPartition[2];
		//cout << "-2" << endl;
		
		if( fieldType == 0 )
		{
			// Read data and create ITL field
			scalarFieldData = ITL_ioutil<SCALAR>::readFieldBinarySerial( inputFieldFile, nDim, (int*)dataSize );
			#ifdef BYTE_SWAP
			int nDataElem = dataSize[0]*dataSize[1]*dataSize[2];
			swap((char *)scalarFieldData, nDataElem, sizeof( SCALAR ));
			#endif
			scalarField = new ITL_field_regular<SCALAR>( scalarFieldData, nDim, lowF, highF );
			//cout << "-1" << endl;

			// Alocate memory for array of partitions (sub-fields)
			//cout << "0: " << nblocks << endl;
			subScalarFieldArray = new ITL_field_regular<SCALAR>[nblocks];

			// Partition field
			//cout << "1" << endl;
			scalarField->partitionField( nPartition, &subScalarFieldArray );
		}
		else if( fieldType == 1 )
		{
			// Read data and create ITL field
			vectorFieldData = ITL_ioutil<VECTOR3>::readFieldBinarySerial( inputFieldFile, nDim, (int*)dataSize );
			#ifdef BYTE_SWAP
			int nDataElem = dataSize[0]*dataSize[1]*dataSize[2];
			swap((char *)vectorFieldData, nDataElem*3, sizeof( SCALAR ));
			#endif
			vectorField = new ITL_field_regular<VECTOR3>( vectorFieldData, nDim, lowF, highF );

			// Alocate memory for array of partitions (sub-fields)
			subVectorFieldArray = new ITL_field_regular<VECTOR3>[nblocks];

			// Partition field
			vectorField->partitionField( nPartition, &subVectorFieldArray );
		}

		// Go through each block and compute entropy
		//cout << "2" << endl;
		globalEntropyList = new float[nblocks];
		
		// Read data for local blocks using MPI-IO
		for( int i=0; i<nblocks; i++ )
		{	
			if( fieldType == 0 )
			{
				// Create bin field (Specify range if required)
				//cout << "3" << endl;
				if( histogramLowEnd != histogramHighEnd )
					histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );
				//cout << "4" << endl;
				histMapper_scalar->computeHistogramBinField_Scalar( (subScalarFieldArray+i), &binField, nBin );

				// Initialize class that can compute entropy
				//cout << "5" << endl;
				globalEntropyComputer_scalar = new ITL_globalentropy<SCALAR>( binField, histogram, nBin );

				// Compute entropy
				//cout << "6" << endl;
				globalEntropyComputer_scalar->computeGlobalEntropyOfField( false );

				// Print global entropy
				//cout << "7" << endl;
				globalEntropyList[i] = globalEntropyComputer_scalar->getGlobalEntropy();

				int lowInt[3], highInt[3];
				subScalarFieldArray[i].getBounds( lowInt, highInt );
				printf( "Block Limit: %d, %d, %d, %d, %d, %d, Global Entropy: %f\n", 
						lowInt[0], highInt[0],
						lowInt[1], highInt[1],
						lowInt[2], highInt[2],
						globalEntropyList[i] );
	
				//cout << "8" << endl;
				if( numProcs == 1 && verboseMode == 1)			
				fprintf( dumpFile, "%d, %d, %d, %d, %d, %d, %f\n", 
						lowInt[0], highInt[0],
						lowInt[1], highInt[1],
						lowInt[2], highInt[2],
						 globalEntropyList[i] );
 

				// clear up memory
				//cout << "9" << endl;
				delete globalEntropyComputer_scalar;
				delete binField;
				binField = NULL;
				//delete [] scalarFieldArray;

				//cout << "10" << endl;
			}
			else if( fieldType == 1 )
			{
				// Create bin field (Specify range if required)
				histMapper_vector->computeHistogramBinField_Vector( (subVectorFieldArray+i), &binField, nBin );

				// Initialize class that can compute entropy
				globalEntropyComputer_vector = new ITL_globalentropy<VECTOR3>( binField, histogram, nBin );

				// Compute entropy
				globalEntropyComputer_vector->computeGlobalEntropyOfField( false );

				// Print global entropy
				int lowInt[3], highInt[3];
				subVectorFieldArray[i].getBounds( lowInt, highInt );
				globalEntropyList[i] = globalEntropyComputer_vector->getGlobalEntropy();
				printf( "Block Limit: %d, %d, %d, %d, %d, %d, Global Entropy: %f\n", 
						lowInt[0], highInt[0],
						lowInt[1], highInt[1],
						lowInt[2], highInt[2],
						globalEntropyList[i] );

				if( numProcs == 1 && verboseMode == 1 )	
					fprintf( dumpFile, "%d, %d, %d, %d, %d, %d, %f\n", 
							lowInt[0], highInt[0],
							lowInt[1], highInt[1],
							lowInt[2], highInt[2],
							 globalEntropyList[i] );
 
				// clear up memory
				delete globalEntropyComputer_vector;
				delete binField;
				binField = NULL;

			}

		}// end for
		
		// Clear up memory
		if( subScalarFieldArray != NULL )delete [] subScalarFieldArray;
		if( scalarField != NULL ) delete scalarField;
		if( subVectorFieldArray != NULL )delete [] subVectorFieldArray;
		if( vectorField != NULL ) delete vectorField;

	}// end if
	else if( parallelOpMode = WITH_DIY )
	{
		// Initialize DIY after initializing MPI
		DIY_Init( nDim, ROUND_ROBIN_ORDER, tot_blocks, &nblocks, dataSize, num_threads, MPI_COMM_WORLD );
		if( verboseMode == 1 )	printf( "Process %d: Number of blocks: %d\n", rank, nblocks );

		// Create the blocking and default assignment
  		// Note in the blocking call that we are not adding extra ghost cells, but we
  		// are sharing boundaries between blocks (share_face = 1)
  		int given[3] = {0, 0, 0};

		// Decompose domain 
  		DIY_Decompose( 1, 0, 0, given );

  		// Allocate memory for pointers that will hold block data
  		MPI_Datatype complex;
  		if( fieldType == 1 )
  		{
  		 	MPI_Type_contiguous( 3,MPI_FLOAT,&complex );
  	  		MPI_Type_commit( &complex );
  		}

		// Allocate memory for pointers that will hold block data
		SCALAR* data[nblocks];
		VECTOR3* vectordata[nblocks];		
		if( fieldType == 0 ) memset( data, 0, sizeof(SCALAR*) * nblocks );
		if( fieldType == 1 ) memset( vectordata, 0, sizeof(VECTOR3*) * nblocks );

		// Serially visit the blocks (?)
		int* diy_min = new int[3*nblocks];
		int* diy_max = new int[3*nblocks];
		int* diy_size = new int[3*nblocks];
       	
 		for (int i = 0; i < nblocks; i++)
		{ 
    		// Allocate memory for block
 			DIY_Block_starts_sizes( i, &diy_min[3*i], &diy_size[3*i] );
 			//if( fieldType == 0 )
 			//	data[i] = new SCALAR[diy_size[3*i] * diy_size[3*i+1] * diy_size[3*i+2]];
 			//else if( fieldType == 1 )
 			//	vectordata[i] = new VECTOR3[diy_size[3*i] * diy_size[3*i+1] * diy_size[3*i+2]];

  			// post a read for the block
    		if( fieldType == 0 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, DIY_FLOAT, (void**)&(data[i]));
     		if( fieldType == 1 ) DIY_Add_data_raw( &diy_min[3*i], &diy_size[3*i], inputFieldFile, complex, (void**)&(vectordata[i]));

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

			//SCALAR m = ITL_util<SCALAR>::Min( data[i], 64*64*64 );
			//SCALAR M = ITL_util<SCALAR>::Max( data[i], 64*64*64 );
			//cout << m << " " << M << endl;
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

		// Allocate memory for storing blockwise global entropies 		
		globalEntropyList = new float[nblocks];

		starttime = ITL_util<float>::startTimer();		
		// Scan through the blocks, compute global entropy and list values 
		int nBlockElem = 0;
		for( int k=0; k<nblocks; k++ )
		{

			lowF[0] = diy_min[3*k];		highF[0] = diy_max[3*k];
			lowF[1] = diy_min[3*k+1];	highF[1] = diy_max[3*k+1];
			lowF[2] = diy_min[3*k+2];	highF[2] = diy_max[3*k+2];

			if( fieldType == 0 )
			{
				// Initialize ITL scalar field with block data
				scalarField = new ITL_field_regular<SCALAR>( data[k], nDim, lowF, highF );
				#ifdef DEBUG_MODE
				SCALAR m = ITL_util<SCALAR>::Min( data[k], 64*64*64 );
				SCALAR M = ITL_util<SCALAR>::Max( data[k], 64*64*64 );
				printf( "Block value range: %g %g\n", m, M );
				#endif
		
				// Create bin field
				histMapper_scalar->computeHistogramBinField_Scalar( scalarField, &binField, nBin );

				// Initialize class that can compute entropy
				globalEntropyComputer_scalar = new ITL_globalentropy<SCALAR>( binField, histogram, nBin );

				// Compute histogram
				globalEntropyComputer_scalar->computeHistogramFrequencies();

				// Compute entropy 
				//globalEntropyComputer_scalar->computeGlobalEntropyOfField( false );

				// Print global entropy
				//globalEntropyList[k] = globalEntropyComputer_scalar->getGlobalEntropy();

				if( verboseMode == 1 ) 
					printf( "Block Limits: %d, %d, %d, %d, %d, %d, Global entropy: %g\n", diy_min[3*k], diy_max[3*k], \
													     diy_min[3*k+1], diy_max[3*k+1], \
													     diy_min[3*k+2], diy_max[3*k+2], globalEntropyList[k] );
				if( numProcs == 1 && verboseMode == 1 )
					fprintf( dumpFile, "%d, %d, %d, %d, %d, %d, %g\n", diy_min[3*k], diy_max[3*k], \
											 diy_min[3*k+1], diy_max[3*k+1], \
											 diy_min[3*k+2], diy_max[3*k+2], globalEntropyList[k] );

				// Clear up
				delete globalEntropyComputer_scalar;
				delete binField;
				binField = NULL;
				delete scalarField;

			}
			if( fieldType == 1 )
			{
				// Initialize ITL vector field with block data
				//cout << "1" << endl;
				vectorField = new ITL_field_regular<VECTOR3>( vectordata[k], nDim, lowF, highF );

				// Create bin field
				//cout << "2" << endl;
				histMapper_vector->computeHistogramBinField_Vector( vectorField, &binField, nBin );

				// Initialize class that can compute entropy
				//cout << "3" << endl;
				globalEntropyComputer_vector = new ITL_globalentropy<VECTOR3>( binField, histogram, nBin );

				// Compute histogram
				globalEntropyComputer_vector->computeHistogramFrequencies();
	
				// Compute entropy
				//cout << "4" << endl;
				//globalEntropyComputer_vector->computeGlobalEntropyOfField( false );

				// Print global entropy
				//cout << "5" << endl;
				//globalEntropyList[k] = globalEntropyComputer_vector->getGlobalEntropy();
				if( verboseMode == 1 ) 
					printf( "Block Limits: %d, %d, %d, %d, %d, %d, Global entropy: %f\n", diy_min[3*k], diy_max[3*k], \
													     diy_min[3*k+1], diy_max[3*k+1], \
													     diy_min[3*k+2], diy_max[3*k+2], globalEntropyList[k] );
				if( numProcs == 1 && verboseMode == 1 )
					fprintf( dumpFile, "%d, %d, %d, %d, %d, %d, %f\n", diy_min[3*k], diy_max[3*k], \
											 diy_min[3*k+1], diy_max[3*k+1], \
											 diy_min[3*k+2], diy_max[3*k+2], globalEntropyList[k] );

				// Clear up
				delete globalEntropyComputer_vector;
				delete binField;
				binField = NULL;
				delete vectorField;
			}

		}// End for loop
		execTime[1] = ITL_util<float>::endTimer( starttime );
	
		// Write global entropy 
		if( verboseMode == 1 ) printf( "Writing blockwise global entropy ...\n" );
		/*
 		//MPI_Datatype *dtype = new MPI_Datatype; // datatype for output
  		//MPI_Type_contiguous( 1, MPI_FLOAT, dtype );
  		//MPI_Type_contiguous( nblocks, MPI_FLOAT, dtype );
		//MPI_Type_commit( dtype );

		MPI_Barrier( MPI_COMM_WORLD ); // everyone synchronizes again
  			
		float **listptr = new float*[nblocks];
		for( int i=0; i<nblocks; i++ )
		{
			listptr[i] = new float[1];
			listptr[i][0] = globalEntropyList[i];
		}
		
		starttime = ITL_util<float>::startTimer();	
		DIY_Write_open_all( outFile, 0 );
  		DIY_Write_blocks_all( (void **)listptr, nblocks, NULL, 0, &CreateType );
  		DIY_Write_close_all();
		execTime[2] = ITL_util<float>::endTimer( starttime );		
		*/

		// Runtime
		if( verboseMode == 1 ) 	printf( "%d: Read/Computation/Write Time: %g, %g, %g seconds\n", rank, execTime[0], execTime[1], execTime[2] );
		else 					printf( "%d, %g, %g, %g\n", rank, execTime[0], execTime[1], execTime[2] );
		
		// Clear up
		delete [] diy_min;
		delete [] diy_max;
		delete [] diy_size;
	
		// Finalize BIL
		if( fieldType == 1 ) MPI_Type_free( &complex );
		DIY_Finalize();

	}// End if-else

	// Close file
	if( dumpFile != NULL )	fclose( dumpFile );

	// Finalize MPI
	MPI_Finalize();

}// End main




