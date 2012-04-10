#include "ITL_entropycore.h"

float
ITL_entropycore::computeEntropy_HistogramBased( int* freqArray, int nPoint, int nBin, bool toNormalize )
{
	// Initialize probability array
	#if defined( _WIN32 ) || defined( _WIN64 )
		float* probArray = new float[nBin];
	#else
		float probArray[nBin];
	#endif
		
	// Turn count into probabilities
	for( int i=0; i<nBin; i++ )
		probArray[i] = freqArray[i] / (float)nPoint;

	// Compute negative entropy
	float entropy = 0;
	for( int i = 0; i<nBin; i++ )
	{
		entropy += ( probArray[i] * ( probArray[i] == 0 ? 0 : ( log( probArray[i] ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;
	
	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (float)nBin ) / log( 2.0f ) );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] probArray;
	#endif

	return entropy;
	
}// End function

float
ITL_entropycore::computeEntropy_HistogramBased2( float* probArray, int nBin, bool toNormalize )
{

	// Compute negative entropy
	float entropy = 0;
	for( int i = 0; i<nBin; i++ )
	{
		entropy += ( probArray[i] * ( probArray[i] == 0 ? 0 : ( log( probArray[i] ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;

	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (float)nBin ) / log( 2.0f ) );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] probArray;
	#endif

	return entropy;

}// End function

double
ITL_entropycore::computeEntropy_HistogramBased2( double* probArray, int nBin, bool toNormalize )
{

	// Compute negative entropy
	double entropy = 0;
	for( int i = 0; i<nBin; i++ )
	{
		entropy += ( probArray[i] * ( probArray[i] == 0 ? 0 : ( log( probArray[i] ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;

	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (double)nBin ) / log( 2.0 ) );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] probArray;
	#endif

	return entropy;

}// End function

float
ITL_entropycore::computeEntropy_HistogramBased( int* binIds, int* freqArray, int nPoint, int nBin, bool toNormalize )
{
	// Initialize probability array
	#if defined( _WIN32 ) || defined( _WIN64 )
		float* probArray = new float[nBin];
	#else
		float probArray[nBin];
	#endif

	for( int i=0; i<nBin; i++ )
	{
		freqArray[i] = 0;
		probArray[i] = 0;
	}

	// Scan through bin Ids and keep count
	for( int i=0; i<nPoint; i++ )
		freqArray[ binIds[i] ] ++;

	// Turn count into probabilities
	for( int i=0; i<nBin; i++ )
		probArray[i] = freqArray[i] / (float)nPoint;

	// Compute negative entropy
	float entropy = 0;
	for( int i = 0; i<nBin; i++ )
	{
		entropy += ( probArray[i] * ( probArray[i] == 0 ? 0 : ( log( probArray[i] ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;
	
	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (float)nBin ) / log( 2.0f ) );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] probArray;
	#endif

	return entropy;
	
}// End function

float ITL_entropycore::computeEntropy_HistogramBased( float* freqArray, int nBin, bool toNormalize )
{
	float total = 0;
	for( int i=0; i<nBin; i++ )
	{
		total += freqArray[i];
	}

	// Turn into probabilities
	for( int i=0; i<nBin; i++ )
		freqArray[i] = freqArray[i] / total;

	// Compute negative entropy
	float entropy = 0;
	for( int i = 0; i<nBin; i++ )
	{
		entropy += ( freqArray[i] * ( freqArray[i] == 0 ? 0 : ( log( freqArray[i] ) / log(2.0) ) ) );
	}

	// Change sign
	entropy = -entropy;

	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (float)nBin ) / log( 2.0f ) );

	return entropy;

}// End function

float ITL_entropycore::computeEntropy_KDEBased( float* data, int nPoint, float h, bool toNormalize )
{
	// Compute mean and variance of data
	float mu = ITL_statutil<float>::Mean( data, nPoint );
	float var = ITL_statutil<float>::Variance( data, nPoint, mu );
	printf( "Mean and variance of the field: %f %f\n", mu, var );
	
	// Initialize probability array
	#if defined( _WIN32 ) || defined( _WIN64 )
		float* probArray = new float[nPoint];
	#else
		float probArray[nPoint];
	#endif
	for( int i=0; i<nPoint; i++ )
		probArray[i] = 0;
	
	// Estimale kernel bandwidth
	if( h == 0 )
		h = 1.06 * var * pow( nPoint,-0.2f ) ;
        printf( "Number of points in the field: %d\n", nPoint );	
	printf( "Kernel bandwidth: %f\n", h );

	// Estimate probabilities at each point
	float sumOfKernels = 0;
	for( int i=0; i<nPoint; i++ )
	{
		sumOfKernels = 0;
		for( int j=0; j<nPoint; j++ )
		{
			sumOfKernels += ITL_entropycore::evaluateKernel( ( data[i] - data[j] ) / h, mu, var ); 
		}
		probArray[i] = sumOfKernels / (nPoint*h);
	}
	
	// Compute negative entropy
	float entropy = 0;
	for( int i = 0; i<nPoint; i++ )
		entropy += ( probArray[i] * ( probArray[i] == 0 ? 0 : ( log( probArray[i] ) / log(2.0) ) ) );
	
	// Change sign
	entropy = -entropy;
	
	// Normalize, if required
	if( toNormalize )
		entropy /= ( log( (float)nPoint ) / log( 2.0f ) );

	#if defined( _WIN32 ) || defined( _WIN64 )
		delete [] probArray;
	#endif

	return entropy;
	
}// End function

float ITL_entropycore::computeEntropy_KDEBased( VECTOR3* data, int nPoint, float h, bool toNormalize )
{
	return 0.0f;
}// End function

float ITL_entropycore::evaluateKernel( float x, float mu, float var )
{
	float exponent = - ( x - mu )/ (2*var);
	return ( 1.0f / sqrt( 2*pi*var ) ) * exp( exponent );	
	
}// End function
