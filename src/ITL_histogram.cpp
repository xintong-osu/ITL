/**
 * @file ITL_histogram.cpp
 * Source file for ITL_histogram.
 * Created on: Nov 19, 2010
 * @author Teng-Yok
 * @author Abon
 */
#include "ITL_histogram.h"

ITL_histogram::ITL_histogram( const char* patchFileName, int nBins )
{
	ITL_init_histogram( patchFileName, nBins );
}

void
ITL_histogram::ITL_init_histogram( const char* patchFileName, int nBins )
{
	// Initialize histogram parameters
	if( patchFileName[0] == '!' )
	{
		printf( "Reading default patch file...\n" );
		ITL_histogram::readPatches_header();
	}
	else
	{
		printf( "Reading patch file ..\n" );		
		ITL_histogram::readPatches_region( patchFileName );
	}	
	// Compute the histogram look up table - once for all
	ITL_histogram::ITL_compute_table( nBins );

}// end function

// Read the lookup table provided in the header file
void
ITL_histogram::readPatches_header()
{
	static char szPatch[] = { 
	#include "ITL_patch.h" 
	};
	
	char *szToken = strtok(szPatch, "\n");
	float f2Temp[2];

	for(int i=0;i<360; i++)
	{
		float fAngle_radian;
		sscanf( szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
		ITL_histogram::theta[i*2+0] = f2Temp[0];
		ITL_histogram::theta[i*2+1] = f2Temp[1];
		szToken = strtok(NULL, "\n");
		//printf( "%f %f\n", f2Temp[0], f2Temp[1] );

		sscanf(szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
		ITL_histogram::phi[i*2+0] = f2Temp[0];
		ITL_histogram::phi[i*2+1] = f2Temp[1];
		szToken = strtok(NULL, "\n");
		//printf( "%f %f\n", f2Temp[0], f2Temp[1] );
	}
}

// read the lookup table that map the thetas and phis to the patch index;
// the lookup table is stored in the file 'patch.txt,' which contains 360 patches.
void
ITL_histogram::readPatches_region( const char* patchFileName )
{
	if( true == ITL_histogram::bIsPatchRead )
		return;

	float x[2],y[2];
	FILE* fp = fopen( patchFileName, "rb" );
	if( fp == NULL )
		printf( "file open error\n" );

	for(int i=0;i<360; i++)
	{
		fscanf(fp,"%f, %f", &x[0],&y[0]);
		fscanf(fp,"%f, %f", &x[1],&y[1]);

		ITL_histogram::theta[i*2+0]=x[0];
		ITL_histogram::theta[i*2+1]=y[0];
		ITL_histogram::phi[i*2+0]=x[1];
		ITL_histogram::phi[i*2+1]=y[1];
	}
	fclose(fp);

	ITL_histogram::bIsPatchRead = true;
}

void ITL_histogram::ITL_compute_table( int nBins )
{
	for( int t = 0; t < ITL_histogram::iNrOfThetas; t++ )
		for( int p = 0; p < ITL_histogram::iNrOfPhis; p++ )
		{
			float fTheta =	M_PI * 2.0f * float(t) / ITL_histogram::fNrOfThetas;
			float fPhi =	M_PI * float(p) / ITL_histogram::fNrOfPhis;
			int iBin = ITL_histogram::get_bin_by_angle(fTheta, fPhi, nBins );
			if( iBin >= 0 )
				ITL_histogram::piAngleMap[t*ITL_histogram::iNrOfPhis + p] = iBin;

		}

	// set the flag to true
	ITL_histogram::bIsAngleMapInitialized = true;

}// end function

// convert the angles in the spherical coordinates (mytheta, myphi) to the patch index
// according to the lookup table composed of theta and phi.
// The #entries in the lookup table is specified by the parameter 'binnum'.
int
ITL_histogram::get_bin_by_angle(float mytheta, float myphi, int binnum )
{
	for(int i=0; i<binnum;i++)
	{
		if( (mytheta>=theta[i*2+0]) && (mytheta<=theta[i*2+1])&&
			(myphi>=phi[i*2+0]) && (myphi<=phi[i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta+2*pi)>=theta[i*2+0]) && ((mytheta+2*pi)<=theta[i*2+1])&&
			(myphi>=phi[i*2+0]) && (myphi<=phi[i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta)>=theta[i*2+0]) && ((mytheta)<=theta[i*2+1])&&
			((myphi+2*pi)>=phi[i*2+0]) && ((myphi+2*pi)<=phi[i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta+2*pi)>=theta[i*2+0]) && ((mytheta+2*pi)<=theta[i*2+1])&&
			((myphi+2*pi)>=phi[i*2+0]) && ((myphi+2*pi)<=phi[i*2+1])
			)
		{
			return i;
		}
	}
	return -1;
}
// ADD-BY-LEETEN 02/02/2010-END

// compute the theta
float
ITL_histogram::getAngle(float x, float y)
{

	if((x==0)&&(y==0))
		return 0;
	else
	{
		return (pi+(atan2(y,x)));
	}
}

// compute phi
float
ITL_histogram::getAngle2(float x, float y)
{

	if((x==0)&&(y==0))
		return 0;
	else
	{
		return fabs(pi/2.0-(atan2(y,x)));
	}
}

int
ITL_histogram::get_bin_number_3D( SCALAR v, int nBins )
{
	return 0;
}

// convert the vector from Cartesian coodinates to the patch index via the specified lookup table
int
ITL_histogram::get_bin_number_3D( VECTOR3 v, int nBins )
{
	// convert the vector from Cartesian coodinates to sphercial coordinates;
	// the magnitude is ignored.
	float mytheta=getAngle(v.x(), v.y());//0~2pi
	float myphi=  getAngle2(sqrt(v.x()*v.x()+v.y()*v.y()), v.z());//0~pi

	if(!ITL_histogram::theta || !ITL_histogram::phi)
	{
		printf("read in the regions first\n");
		return -1;
	}

	if( ITL_histogram::bIsAngleMapInitialized == false )
		ITL_histogram::ITL_compute_table( nBins );

	// map the angle to the bin number
	int iTheta	=	min(ITL_histogram::iNrOfThetas - 1,	int(ITL_histogram::fNrOfThetas * mytheta	/ (M_PI * 2.0)));
	int iPhi	=	min(ITL_histogram::iNrOfPhis - 1,		int(ITL_histogram::fNrOfPhis * myphi	/ M_PI));

	//cout << ITL_histogram::piAngleMap[ iTheta * ITL_histogram::iNrOfPhis + iPhi] << endl;
	return ITL_histogram::piAngleMap[ iTheta * ITL_histogram::iNrOfPhis + iPhi];
}
