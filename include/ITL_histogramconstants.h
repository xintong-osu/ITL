/**
 * Histogram parameters.
 * Created on: Nov 21, 2010
 * @author Teng-Yok
 * @author Abon
 * @see ITL_histogram
 */

#ifndef ITL_HISTOGRAMCONSTANTS_H_
#define ITL_HISTOGRAMCONSTANTS_H_

#include "ITL_histogram.h"

#if	0	// DEL-BY-LEETEN 07/18/2011-BEGIN
	// Initialize static parameters of a histogram
	bool ITL_histogram::bIsAngleMapInitialized = false;
	bool ITL_histogram::bIsPatchRead = false;
	int ITL_histogram::iNrOfThetas = 720;
	int ITL_histogram::iNrOfPhis =	360;
	float ITL_histogram::fNrOfThetas = float(ITL_histogram::iNrOfThetas);
	float ITL_histogram::fNrOfPhis	= float(ITL_histogram::iNrOfPhis);
	int* ITL_histogram::piAngleMap = new int[ITL_histogram::iNrOfThetas*ITL_histogram::iNrOfPhis];
	float* ITL_histogram::theta=new float[2*360];
	float* ITL_histogram::phi=new float[2*360];
#endif	// DEL-BY-LEETEN 07/18/2011-END

#endif /* ITL_HISROGRAMCONSTANTS_H_ */
