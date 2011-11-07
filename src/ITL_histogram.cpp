/**
 * @file ITL_histogram.cpp
 * Source file for ITL_histogram.
 * Created on: Nov 19, 2010
 * @author Teng-Yok
 * @author Abon
 */
#include "ITL_histogram.h"

// ADD-BY-LEETEN 07/18/2011-BEGIN
// Initialize static parameters of a histogram
//bool ITL_histogram::bIsAngleMapInitialized = false;
//bool ITL_histogram::bIsPatchRead = false;
//int ITL_histogram::iNrOfThetas = 720;
//int ITL_histogram::iNrOfPhis =	360;
//float ITL_histogram::fNrOfThetas = float(ITL_histogram::iNrOfThetas);
//float ITL_histogram::fNrOfPhis	= float(ITL_histogram::iNrOfPhis);
//int* ITL_histogram::piAngleMap = new int[ITL_histogram::iNrOfThetas*ITL_histogram::iNrOfPhis];
//float* ITL_histogram::theta=new float[2*360];
//float* ITL_histogram::phi=new float[2*360];
// ADD-BY-LEETEN 07/18/2011-END


ITL_histogram::ITL_histogram( const char* patchFileName, int nbin )
{
	nResolution = 1;

	// Initialize variables
	nBin[0] = nbin;
	iNrOfThetas[0] = nBin[0]*2;
	iNrOfPhis[0] = nBin[0];
	fNrOfThetas[0] = float(iNrOfThetas[0]);
	fNrOfPhis[0] = float(iNrOfPhis[0]);
	piAngleMap[0] = new int[iNrOfThetas[0]*iNrOfPhis[0]];
	theta[0] = new float[2*nBin[0]];
	phi[0] = new float[2*nBin[0]];
	
	// Initialize bin map
	ITL_init_histogram( patchFileName );
	
	// Set flag
	bIsPatchRead = true;

}

ITL_histogram::ITL_histogram( const char** patchfilename, int *nbin, int nresolution )
{
	assert( nresolution <= 5 );
	nResolution = nresolution;
	
	// Initialize variables
	for( int iR=0; iR<nResolution; iR++ )
	{
		nBin[iR] = nbin[iR];
		iNrOfThetas[iR] = nBin[iR]*2;
		iNrOfPhis[iR] = nBin[iR];
		fNrOfThetas[iR] = float(iNrOfThetas[iR]);
		fNrOfPhis[iR] = float(iNrOfPhis[iR]);
		piAngleMap[iR] = new int[iNrOfThetas[iR]*iNrOfPhis[iR]];	
		theta[iR] = new float[2*nBin[iR]];
		phi[iR] = new float[2*nBin[iR]];
	
		// Initialize bin map
		ITL_init_histogram( patchfilename[iR], iR );
	}
	
	// Set flag
	bIsPatchRead = true;
}

void
ITL_histogram::ITL_init_histogram( const char* patchFileName, int iRes )
{
	// Initialize histogram parameters
	if( patchFileName[0] == '!' )
	{
		// The patch provided in the header is only for
		// one resolution with 360 bins
		assert( iRes == 0 );
		
		//printf( "Reading default patch file...\n" );
		readPatches_header();
	}
	else
	{
		//printf( "Reading patch file ..\n" );		
		readPatches_region( patchFileName, iRes );
	}	
	
	// Compute the histogram look up table - once for all
	ITL_compute_table();

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

	for(int i=0;i<nBin[0]; i++)
	{
		float fAngle_radian;
		sscanf( szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
		theta[0][i*2+0] = f2Temp[0];
		theta[0][i*2+1] = f2Temp[1];
		szToken = strtok(NULL, "\n");
		//printf( "%f %f\n", f2Temp[0], f2Temp[1] );

		sscanf(szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
		phi[0][i*2+0] = f2Temp[0];
		phi[0][i*2+1] = f2Temp[1];
		szToken = strtok(NULL, "\n");
		//printf( "%f %f\n", f2Temp[0], f2Temp[1] );
	}
}

// read the lookup table that map the thetas and phis to the patch index;
// the lookup table is stored in the file 'patch.txt,' which contains 360 patches.
void
ITL_histogram::readPatches_region( const char* patchFileName, int iRes )
{
	//if( true == ITL_histogram::bIsPatchRead )
	//	return;

	float x[2],y[2];
	cout << "Reading patch file" << patchFileName << endl;
	FILE* fp = fopen( patchFileName, "r" );
	if( fp == NULL )
		printf( "file open error\n" );

	for(int i=0;i<nBin[iRes]; i++)
	{
		fscanf(fp,"%f, %f", &x[0],&y[0]);
		fscanf(fp,"%f, %f", &x[1],&y[1]);

		theta[iRes][i*2+0]=x[0];
		theta[iRes][i*2+1]=y[0];
		phi[iRes][i*2+0]=x[1];
		phi[iRes][i*2+1]=y[1];
	}
	fclose(fp);
	cout << "Reading patch file done" << endl;
	//ITL_histogram::bIsPatchRead = true;
}

void
ITL_histogram::ITL_compute_table()
{
	for( int iR=0; iR<nResolution; iR++ )
	{
		ITL_compute_table( iR );
	}
}

void
ITL_histogram::ITL_compute_table( int iRes )
{
	for( int t = 0; t < iNrOfThetas[iRes]; t++ )
		for( int p = 0; p < iNrOfPhis[iRes]; p++ )
		{
			float fTheta =	M_PI * 2.0f * float(t) / fNrOfThetas[iRes];
			float fPhi =	M_PI * float(p) / fNrOfPhis[iRes];
			int iBin = get_bin_by_angle(fTheta, fPhi, nBin[iRes], iRes );
			if( iBin >= 0 )
				piAngleMap[iRes][t*iNrOfPhis[iRes] + p] = iBin;
		}

	// set the flag to true
	bIsAngleMapInitialized[iRes] = true;

}// end function

// convert the angles in the spherical coordinates (mytheta, myphi) to the patch index
// according to the lookup table composed of theta and phi.
// The #entries in the lookup table is specified by the parameter 'binnum'.
int
ITL_histogram::get_bin_by_angle( float mytheta, float myphi, int binnum, int iRes )
{
	for(int i=0; i<binnum;i++)
	{
		if( (mytheta>=theta[iRes][i*2+0]) && (mytheta<=theta[iRes][i*2+1])&&
			(myphi>=phi[iRes][i*2+0]) && (myphi<=phi[iRes][i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta+2*pi)>=theta[iRes][i*2+0]) && ((mytheta+2*pi)<=theta[iRes][i*2+1])&&
			(myphi>=phi[iRes][i*2+0]) && (myphi<=phi[iRes][i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta)>=theta[iRes][i*2+0]) && ((mytheta)<=theta[iRes][i*2+1])&&
			((myphi+2*pi)>=phi[iRes][i*2+0]) && ((myphi+2*pi)<=phi[iRes][i*2+1])
			)
		{
			return i;
		}
	}
	for(int i=0; i<binnum;i++)
	{
		if( ((mytheta+2*pi)>=theta[iRes][i*2+0]) && ((mytheta+2*pi)<=theta[iRes][i*2+1])&&
			((myphi+2*pi)>=phi[iRes][i*2+0]) && ((myphi+2*pi)<=phi[iRes][i*2+1])
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
ITL_histogram::get_bin_number_3D( VECTOR3 v, int iRes )
{
	// convert the vector from Cartesian coodinates to sphercial coordinates;
	// the magnitude is ignored.
	float mytheta=getAngle(v.x(), v.y());//0~2pi
	float myphi=  getAngle2(sqrt(v.x()*v.x()+v.y()*v.y()), v.z());//0~pi

	if( !theta[iRes] || !phi[iRes] )
	{
		printf("read in the regions first\n");
		return -1;
	}

	if( bIsAngleMapInitialized[iRes] == false )
		ITL_compute_table( iRes );

	// map the angle to the bin number
	int iTheta	=	min( iNrOfThetas[iRes] - 1,	int( fNrOfThetas[iRes] * mytheta / (M_PI * 2.0)));
	int iPhi	=	min( iNrOfPhis[iRes] - 1, int( fNrOfPhis[iRes] * myphi	/ M_PI));

	return piAngleMap[iRes][ iTheta * iNrOfPhis[iRes] + iPhi];
}

// convert the 2D vector from Cartesian coodinates to the patch index via the specified lookup table
int
ITL_histogram::get_bin_number_2D( VECTOR3 v, int nbin )
{
	// Convert the vector from Cartesian coodinates to sphercial coordinates;
	// the magnitude is ignored.
	float mytheta= getAngle(v.x(), v.y());//0~2pi

	return (int)floor( (mytheta / (2*pi) ) * nbin );
}

