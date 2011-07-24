/*
 * CSphereSpace.cpp
 *
 *  Created on: Jul 21, 2011
 *      Author: leeten
 */

#include <mpi.h>
#include "liblog.h"
#include "ITL_SphereSpace.h"

bool CSphereSpace::bIsInitialized = false;
TBuffer<CRange> CSphereSpace::pcDefaultThetaBins;
TBuffer<CRange> CSphereSpace::pcDefaultPhiBins;

//! load the default mapping
void
CSphereSpace::_LoadDefaultMapping()
{
	if( true == bIsInitialized )
		return;

	static char szDefaultMapping[] = {
		#include "ITL_patch.h"
	};

	pcDefaultThetaBins.New(DEFAULT_NR_OF_PATCHES);
	pcDefaultPhiBins.New(DEFAULT_NR_OF_PATCHES);

	char *szToken = strtok(szDefaultMapping, "\n");
	for(int i = 0; i < DEFAULT_NR_OF_PATCHES; i++)
	{
		float f2Theta[2];
		sscanf( szToken,"%f,%f", &f2Theta[0], &f2Theta[1] );
		pcDefaultThetaBins[i]._Set((double)f2Theta[0], (double)f2Theta[1] );
		szToken = strtok(NULL, "\n");

		float f2Phi[2];
		sscanf(szToken,"%f,%f", &f2Phi[0], &f2Phi[1] );
		pcDefaultPhiBins[i]._Set((double)f2Phi[0], (double)f2Phi[1] );
		szToken = strtok(NULL, "\n");
	}
	bIsInitialized = true;
}

int
// MOD-BY-LEETEN 07/23/2011-FROM:
// CSphereSpace::IGetNrOfPatches()
// TO:
CSphereSpace::IGetNrOfPatches
(
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	return iNrOfPatches;
}

int
CSphereSpace::IGetNrOfThetas
(
// MOD-BY-LEETEN 07/23/2011-FROML
// )
// TO:
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	return this->p2DiMapping.iHeight;
}

int
CSphereSpace::IGetNrOfPhis
(
// MOD-BY-LEETEN 07/23/2011-FROM:
// )
// TO:
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	return this->p2DiMapping.iWidth;
}

// compute the theta
double
// MOD-BY-LEETEN 07/23/2011-FROM:
// CSphereSpace::DGetAngle(double x, double y)
// TO:
CSphereSpace::DGetAngle
(
 double x, 
 double y
) const
// MOD-BY-LEETEN 07/23/2011-END
{

	if( 0.0 == x && 0.0 == y )
		return 0.0;
	else
		return M_PI + atan2(y, x);
}

// compute phi
double
// MOD-BY-LEETEN 07/23/2011-FROM:
// CSphereSpace::DGetAngle2(double x, double y)
// TO:
CSphereSpace::DGetAngle2
(
 double x, 
 double y
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	if( 0.0 == x && 0.0 == y )
		return 0.0;
	else
		return fabs(M_PI/2.0 - atan2(y, x));
}


int
// MOD-BY-LEETEN 07/23/2011-FROM:
// CSphereSpace::IMapVectorToPatch(double pdVector[])
// TO:
CSphereSpace::IMapVectorToPatch
(
 const double pdVector[]
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	// convert the vector from Cartesian coodinates to sphercial coordinates;
	// the magnitude is ignored.
	double mytheta = DGetAngle(pdVector[0], pdVector[1]);//0~2pi
	double myphi = 	 DGetAngle2(sqrt(pdVector[0] * pdVector[0] + pdVector[1] * pdVector[1]), pdVector[2]);//0~pi

	// map the angle to the bin number
	int iNrOfThetas = IGetNrOfThetas();
	int iNrOfPhis = IGetNrOfPhis();
	int iTheta	=	max(0, min(iNrOfThetas 	- 1,	(int)floor(mytheta *(double)iNrOfThetas / (M_PI * 2.0))));
	int iPhi	=	max(0, min(iNrOfPhis 	- 1,	(int)floor(myphi * 	(double)iNrOfPhis 	/ M_PI)));

	//cout << ITL_histogram::piAngleMap[ iTheta * ITL_histogram::iNrOfPhis + iPhi] << endl;
	return this->p2DiMapping.GetElement(iPhi, iTheta);
}

// convert the angles in the spherical coordinates (mytheta, myphi) to the patch index
// according to the lookup table composed of theta and phi.
// The #entries in the lookup table is specified by the parameter 'binnum'.
int
// MOD-BY-LEETEN 07/23/2011-FROM:
// CSphereSpace::IGetBinByAngle(double mytheta, double myphi)
// TO:
CSphereSpace::IGetBinByAngle
(
     double mytheta, 
     double myphi
) const
// MOD-BY-LEETEN 07/23/2011-END
{
	for(int i = 0; i < iNrOfPatches; i++)
		if(
			( 	pcThetaBins[i].dMin	<= 	mytheta && 	mytheta <= pcThetaBins[i].dMax &&
				pcPhiBins[i].dMin 	<= 	myphi  && 	myphi 	<= pcPhiBins[i].dMax ) ||
			( 	pcThetaBins[i].dMin <= 	mytheta + M_PI * 2 && 	mytheta + M_PI * 2 <= pcThetaBins[i].dMax &&
				pcPhiBins[i].dMin <= 	myphi  && 	myphi <= pcPhiBins[i].dMax ) ||
			( 	pcThetaBins[i].dMin <= 	mytheta && 	mytheta <= pcThetaBins[i].dMax &&
				pcPhiBins[i].dMin <= 	myphi  	+ M_PI * 2 && 	myphi 	+ M_PI * 2 <= pcPhiBins[i].dMax ) ||
			( 	pcThetaBins[i].dMin <= 	mytheta + M_PI * 2 && 	mytheta + M_PI * 2 <= pcThetaBins[i].dMax &&
				pcPhiBins[i].dMin <= 	myphi   + M_PI * 2 && 	myphi 	+ M_PI * 2 <= pcPhiBins[i].dMax ) )
			return i;

	// Get the rank of the current processors
	int iRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
	if(!iRank)
	{
		for(int i = 0; i < iNrOfPatches; i++)
			LOG_ERROR(fprintf(stderr, "%d: %e - %e, %e - %e.", i,
					pcThetaBins[i].dMin, pcThetaBins[i].dMax,
					pcPhiBins[i].dMin, pcPhiBins[i].dMax) );
		ASSERT_OR_LOG(false, fprintf(stderr, "%e, %e: invalid theta & phi.", mytheta, myphi));
	}
	return -1;
}

//! the method to load the default mapping
void
CSphereSpace::_CopyDefaultMapping()
{
  // ADD-BY-LEETEN 07/23/2011-BEGIN
  _LoadDefaultMapping();
  // ADD-BY-LEETEN 07/23/2011-END
	this->iNrOfPatches = DEFAULT_NR_OF_PATCHES;
	this->pcThetaBins.New(DEFAULT_NR_OF_PATCHES);
	this->pcPhiBins.New(DEFAULT_NR_OF_PATCHES);

	// ADD-BY-LEETEN 07/23/2011-BEGIN
	memcpy(&pcThetaBins[0],	&pcDefaultThetaBins[0],	sizeof(pcThetaBins[0]) * pcThetaBins.USize());
	memcpy(&pcPhiBins[0], 	&pcDefaultPhiBins[0], 	sizeof(pcPhiBins[0])   * pcPhiBins.USize()	);
	// create the table to accelerate the mapping
	_ComputeMapping(DEFAULT_NR_OF_THETA_SAMPLES, DEFAULT_NR_OF_PHI_SAMPLES);
	// ADD-BY-LEETEN 07/23/2011-END
}

void
CSphereSpace::_ComputeMapping
(
	int iNrOfThetas,
	int iNrOfPhis
)
{
	p2DiMapping.alloc(iNrOfPhis, iNrOfThetas);
	for(int i = 0, 	t = 0; t < iNrOfThetas; t++ )
		for( int 	p = 0; p < iNrOfPhis; 	p++, i++ )
		{
			float fTheta =	float(M_PI * 2.0 * double(t) / (double)iNrOfThetas);
			float fPhi =	float(M_PI * double(p) / (double)iNrOfPhis);
			int iBin = IGetBinByAngle(fTheta, fPhi);
			if( iBin >= 0 )
				p2DiMapping[i] = iBin;
		}
}

CSphereSpace::CSphereSpace() {
	// TODO Auto-generated constructor stub
	CSphereSpace::_LoadDefaultMapping();
}

CSphereSpace::~CSphereSpace() {
	// TODO Auto-generated destructor stub
}
