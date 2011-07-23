/*
 * CSphereSpace.h
 *
 *  Created on: Jul 21, 2011
 *      Author: leeten
 */

#ifndef CSPHERESPACE_H_
#define CSPHERESPACE_H_

#include <vector>
using namespace std;
#include <math.h>

#include "ITL_Range.h"

#include "liblog.h"
#include "libbuf.h"
#include "libbuf2d.h"
#include "libbuf3d.h"

class CSphereSpace
{
	static bool bIsInitialized;
	static TBuffer<CRange> pcDefaultThetaBins;
	static TBuffer<CRange> pcDefaultPhiBins;

	int iNrOfPatches;
	TBuffer<CRange> pcThetaBins;
	TBuffer<CRange> pcPhiBins;

	TBuffer2D<int> p2DiMapping;
public:
	enum {
		DEFAULT_NR_OF_THETA_SAMPLES	= 512,
		DEFAULT_NR_OF_PHI_SAMPLES	= 256,
		DEFAULT_NR_OF_PATCHES		= 360
	};

	int IGetNrOfThetas();
	int IGetNrOfPhis();
	int IGetNrOfPatches();

	void
	_ComputeMapping
	(
		int iNrOfThetas,
		int iNrOfPhis
	);

	int
	IMapVectorToPatch
	(
		double pdVector[]
	);

	// compute the theta
	double
	DGetAngle
	(
		double x,
		double y
	);

	// compute phi
	double
	DGetAngle2
	(
		double x,
		double y
	);

	// convert the angles in the spherical coordinates (mytheta, myphi) to the patch index
	// according to the lookup table composed of theta and phi.
	// The #entries in the lookup table is specified by the parameter 'binnum'.
	int
	IGetBinByAngle
	(
		double mytheta,
		double myphi
	);

	//! the method to load the default mapping
	void _CopyDefaultMapping
	(
	);

	static
	void _LoadDefaultMapping();

	CSphereSpace();
	virtual ~CSphereSpace();
};

#endif /* CSPHERESPACE_H_ */
