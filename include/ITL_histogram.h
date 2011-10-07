/**
 * Histogram class.
 * The histogram assigns a field value to a bin of the histogram.
 * Created on: Nov 17, 2010.
 * @author Teng-Yok
 * @author Abon
 * @see ITL_histogramconstants
 */

#ifndef ITL_HISTOGRAM_H_
#define ITL_HISTOGRAM_H_

#include "ITL_header.h"
#include "ITL_vectormatrix.h"

class ITL_histogram
{
public:
	// ADD-BY-LEETEN 10/07/2011-BEGIN
	enum {
		DEFAULT_NR_OF_BINS = 360
	};
	// ADD-BY-LEETEN 10/07/2011-END

	static float* theta;				 /**< Angle variable. Angle related to spherical coordinates. */
	static float* phi;					 /**< Angle variable. Angle related to spherical coordinates. */

	static bool bIsAngleMapInitialized;	 /**< Boolean flag. Flag is set if the angle maps have been initialized. */
	static bool bIsPatchRead;			 /**< Boolean flag. Flag is set if the patch file has been read. */
	static int iNrOfThetas;				 /**< Histogram parameter. Number of thetas.  */
	static int iNrOfPhis;				 /**< Histogram parameter. Number of phis.  */
	static float fNrOfThetas;			 /**< Histogram parameter. Number of phis.  */
	static float fNrOfPhis;				 /**< Histogram parameter. Number of phis.  */
	static int* piAngleMap;				 /**< Histogram parameter. Mapping from vector to a region in sperical coordinates.  */


public:

	/**
	 * Constructor.
	 * @param patchFileName Path of patch file on disc.
	 * @param nBins Desired number of bins in the histogram.
	 */
	ITL_histogram( const char* patchFileName, int nBins );

	/**
	 * Histogram initialization function.
	 * @param patchFileName Path of patch file on disc.
	 * @param nBins Desired number of bins in the histogram.
	 */
	static void ITL_init_histogram( const char* patchFileName, int nBins );
	/**
	 * Patch header reader. Reads patch information from the in-built header.
	 * 
	 */
	static void readPatches_header();
	/**
	 * Patch file reader.
	 * @param patchFileName Path of patch file on disc.
	 *
	 */
	static void readPatches_region( const char* patchFileName );
	/**
	 * Actual angle map computation doen here.
	 * @param nBins Desired number of bins in the histogram.
	 */
	static void ITL_compute_table( int nBins );

	/**
	 * Actual angle map computation done here.
	 * @param mytheta theta corresponding to the local vector.
	 * @param myphi phi corresponding to the local vector.
	 * @param binnum Desired number of bins in the histogram.
	 */
	static int get_bin_by_angle(float mytheta, float myphi, int binnum );

	/**
	 * Scala-to-bin conversion Routine.
	 * Function converts the scalar to bin id.
	 * @param v at some Cartesian co-ordinate of the field
	 * @param nBins Desired number of bins in the histogram.
	 */
	static int get_bin_number_3D(SCALAR v, int nBins );
	
	/**
	 * Vector-to-bin conversion Routine.
	 * Function converts the vector from Cartesian coodinates to the patch index
	 * via the specified lookup table.
	 * @param v at some Cartesian co-ordinate of the field
	 * @param nBins Desired number of bins in the histogram.
	 */
	static int get_bin_number_3D(VECTOR3 v, int nBins );

	/**
	 * 2D Vector-to-bin conversion Routine.
	 * Function converts the vector from Cartesian coodinates to the patch index
	 * via the specified lookup table.
	 * @param v at some Cartesian co-ordinate of the field
	 * @param nBins Desired number of bins in the histogram.
	 */
	static int get_bin_number_2D(VECTOR3 v, int nBins );


	/**
	 * Theta computation.
	 * @param x x-component of vector
	 * @param y y-computation of vector
	 * @return Theta
	 */
	static float getAngle(float x, float y);

	/**
	 * Phi computation.
	 * @param x x-component of vector
	 * @param y y-computation of vector
	 * @return Phi
	 */
	static float getAngle2(float x, float y);

};

#endif
/* ITL_HISTOGRAM_H_ */
