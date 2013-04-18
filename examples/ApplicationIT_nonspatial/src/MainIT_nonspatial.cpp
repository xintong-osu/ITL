/**
 * @file MainIT_diy_global.cpp
 * Application program for global entropy computation using diy framework
 * Created on: June 20, 2011
 * @author Abon
 * @author Teng-Yok
 */

//#define BYTE_SWAP

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

int* binField = NULL;

ITL_histogram *histogram = NULL;
ITL_histogrammapper<SCALAR>* histMapper_scalar = NULL;
ITL_histogrammapper<VECTOR3>* histMapper_vector = NULL;

int nBin = 100;

double execTime[5];
clock_t starttime, endtime;

/**
 * Main function.
 * Program starts from here.
 */
int
main( int argc, char** argv )
{
	int nDim = 3, nDataPoint = 0;
	int dataSize[3];

	int fieldType = -1;
	int verboseMode = 1;

	float histogramLowEnd = 0.0f;
	float histogramHighEnd = 0.0f;

	float globalEntropy = 0.0f;

	// Read file containing all command line arguments
	ITL_util<float>::getArgs( argv[1], &argNames, &argValues );

	// Parse command line arguments
	inputFieldFile = ITL_util<float>::getArgWithName( "inputField", &argNames, &argValues );
	patchFile =  ITL_util<float>::getArgWithName( "patchFile", &argNames, &argValues );
	fieldType = atoi( ITL_util<float>::getArgWithName( "fieldType", &argNames, &argValues ) );
	nDim = atoi( ITL_util<float>::getArgWithName( "nDim", &argNames, &argValues ) );
	nBin = atoi( ITL_util<float>::getArgWithName( "nBin", &argNames, &argValues ) );
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

	if( fieldType == 0 )
	{
		ITL_ioutil<float>::readFieldBinarySerial( inputFieldFile, nDim, dataSize, &scalarFieldData );

		nDataPoint = dataSize[0]*dataSize[1]*dataSize[2];

		binField = new int[nDataPoint];

		//cout << "-3" << endl;
		if( histogramLowEnd != histogramHighEnd )
			histMapper_scalar->setHistogramRange( histogramLowEnd, histogramHighEnd );

		// Create bin field (Local analysis)
		//cout << "-2" << endl;
		histMapper_scalar->computeHistogramBinField_Scalar_NS( scalarFieldData, nDataPoint, binField, nBin );

		//int bm = ITL_util<int>::Min( binField, nDataPoint );
		//int bM = ITL_util<int>::Max( binField, nDataPoint );
		//cout << bm << " " << bM << endl;

		// Compute frequencies
		//cout << "0" << endl;
		float freqList[nBin];
		memset( freqList, 0.0, sizeof(float)*nBin );
		histMapper_scalar->computeHistogramFrequencies( binField, freqList, nDataPoint );

		// Print / Save in file
		for (int i = 0; i < nBin; i++)
			fprintf(stderr, "freq[%d] = %g\n", i, freqList[i] );
	}
	if( fieldType == 1 )
	{
		ITL_ioutil<VECTOR3>::readFieldBinarySerial( inputFieldFile, nDim, dataSize, &vectorFieldData );

		nDataPoint = dataSize[0]*dataSize[1]*dataSize[2];

		binField = new int[nDataPoint];

		// Create bin field (local analysis)
		//cout << "-3" << endl;
		histMapper_vector->computeHistogramBinField_Vector_NS( vectorFieldData, nDataPoint, binField, nBin );

		// Compute frequencies
		//cout << "0" << endl;
		float freqList[nBin];
		memset( freqList, 0.0, sizeof(float)*nBin );
		histMapper_scalar->computeHistogramFrequencies( binField, freqList, nDataPoint );

		// Print / Save in file
		for (int i = 0; i < nBin; i++)
			fprintf(stderr, "freq[%d] = %g\n", i, freqList[i]);

	}

	// Clear up
	delete [] binField;


}// End main




