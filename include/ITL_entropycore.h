/**
 * Entropy computation core class.
 * Created on: May 03, 2011.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_ENTROPYCORE_H_
#define ITL_ENTROPYCORE_H_

#include "ITL_header.h"
#include "ITL_statutil.h"

class ITL_entropycore
{
public:

public:

	/**
	 * Histogram based entropy computation function.
	 * @param nPoint Number of points. Same as the length of bin array.
	 * @param nBin Number of bins used in histogram computation.
	 */
	static float computeEntropy_HistogramBased( int* freqArray, int nPoint, int nBin, bool toNormalize );

	/**
	 * Histogram based entropy computation function.
	 * @param binIds Pointer to array containing histogram bin assignments of the points.
	 * @param nPoint Number of points. Same as the length of bin array.
	 * @param nBin Number of bins used in histogram computation.
	 */
	static float computeEntropy_HistogramBased( int* binIds, int* freqArray, int nPoint, int nBin, bool toNormalize );

	
	/**
	 * KDE based entropy computation function.
	 * @param data Pointer to data array.
	 * @param nPoint Number of points to be sampled. Equivalent to nBin.
	 * @param h Kernel Bandwidth (Put 0 to allow the software decide the bandwidth)
	 * @param toNormalize TRUE indicates the computed entropy will be normalized.
	 */
	static float computeEntropy_KDEBased( float* data, int nPoint, float h, bool toNormalize );

	/**
	 * KDE based ***Dummy*** entropy computation function.
	 * @param data Pointer to data array.
	 * @param nPoint Number of points to be sampled. Equivalent to nBin.
	 * @param h Kernel Bandwidth (Put 0 to allow the software decide the bandwidth)
	 * @param toNormalize TRUE indicates the computed entropy will be normalized.
	 */
	static float computeEntropy_KDEBased( VECTOR3* data, int nPoint, float h, bool toNormalize );

	/**
	 * Gaussian kernel evaluation function.
	 * Returns N( (x-mu)/sigma |0, 1 ) for some sample x..
	 * @param mu Mean of the point set
	 * @paran var Variance of the point set
	 */
	static float evaluateKernel( float x, float mu, float var );
	
};

#endif
/* ITL_ENTROPYCORE_H_ */
