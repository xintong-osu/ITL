/**
 * @file ITL_histogram.cpp
 * Source file for ITL_histogram.
 * Created on: Nov 19, 2010
 * @author Teng-Yok
 * @author Abon
 */
#include "ITL_histogram.h"

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
	thetaCenter[0] = new float[nBin[0]];
	phiCenter[0] = new float[nBin[0]];
	binCenter[0] = new VECTOR3[nBin[0]];

	// ADD-BY-Tzu-Hsuan-BEGIN-02/13
	binDatas = NULL;
	hist_RangeSet = false;
	// ADD-BY-Tzu-Hsuan-END-02/13
	
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
		thetaCenter[iR] = new float[nBin[iR]];
		phiCenter[iR] = new float[nBin[iR]];
		binCenter[iR] = new VECTOR3[nBin[iR]];

		// ADD-BY-Tzu-Hsuan-BEGIN-02/13
		binDatas = NULL;
		hist_RangeSet = false;
		// ADD-BY-Tzu-Hsuan-END-02/13
	
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

		thetaCenter[0][i] = ( theta[0][i*2+0] + theta[0][i*2+1] ) / 2.0f;
		phiCenter[0][i] = ( phi[0][i*2+0] + phi[0][i*2+1] ) / 2.0f;
		binCenter[0][i] = getXYZ( thetaCenter[0][i], phiCenter[0][i] );
		//binCenter[0][i].Normalize();

	}// end for
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

		thetaCenter[iRes][i] = ( theta[0][i*2+0] + theta[0][i*2+1] ) / 2.0f;
		phiCenter[iRes][i] = ( phi[0][i*2+0] + phi[0][i*2+1] ) / 2.0f;
		binCenter[iRes][i] = getXYZ( thetaCenter[iRes][i], phiCenter[iRes][i] );
		//binCenter[iRes][i].Normalize();

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

// Compute x, y, z
VECTOR3
ITL_histogram::getXYZ( float theta, float phi )
{
	//VECTOR3 v( 1.0f*sin(theta)*cos(phi), 1.0f*sin(theta)*sin(phi), 1.0f*cos(theta) );
	VECTOR3 v( 1.0f*sin(phi)*cos(theta), 1.0f*sin(theta)*sin(phi), 1.0f*cos(phi) );
	return v;
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

	//ADD-BY-ABON-02/26/12-BEGIN
	int binID = (int)floor( (mytheta / (2*pi) ) * nbin );
	#ifdef DEBUG_MODE
	if( binID < 0 || binID >= nbin )
		cout << "binid out of bound ..." << endl;
	#endif
	binID = ITL_util<int>::clamp( binID, 0, nbin-1 );
	//ADD-BY-ABON-02/26/12-END

	return binID;
}

VECTOR3
ITL_histogram::getBinCenter( int binId, int iRes )
{
	return binCenter[iRes][binId];
}

/**
  * Speedup Cross-validation function.
  */

int
ITL_histogram::crossValidateSpeedUp( ITL_field_regular<SCALAR> *scalarField, char *fieldType, int nMax, int start, int step)
{
	// Assume method is 0 (histogram-based)
	// Allocate memory

	int nIter = log((float)(nMax/start))/log(2.0) + 1;
	printf("nIter:%d\n",nIter);
	SCALAR minValue,maxValue;
	double *scoreList = new double[nIter];
	int *nBinArray = new int[nIter];
	int cross_nBin;

	double *freqListPtr = NULL;
	double **freqListPtr2 = NULL;

	for( int i = 0; i<nIter; i++ )
	{
		nBinArray[nIter-i-1] = (int)((float)nMax/powf(2.0f,(float)i));	// MOD-BY-LEETEN 04/09/2012-FROM:		nBinArray[nIter-i-1] = (int)(nMax/(pow(2.0,(float)i)));
	}

	for(int i=0;i<nIter;i++)
		printf("nBinArray[i]:%d\n",nBinArray[i]);

	// Iterate
	double minScore = 10000000;
	int minScoreIndex = -1;
	float h = 0;
	//int N = scalarField->grid->nVertices;
	int N = scalarField->getSize();
	cout << "N " << N << endl;
	//FILE *crossV;
	//crossV = fopen("crossV.txt","w+");
	double firstPart = 0;
	double secondPart = 0;
	double sumOfSquares = 0;

	// Compute the range over which histogram computation needs to be done
	if( hist_RangeSet == false )
	{
		// Get min-max values of the scalar field
		// ******************* Need to resolve typecast issue begin *************************
		//minValue = ITL_util<SCALAR>::Min( (SCALAR*)scalarField->datastore->array, scalarField->grid->nVertices );
		//maxValue = ITL_util<SCALAR>::Max( (SCALAR*)scalarField->datastore->array, scalarField->grid->nVertices );
		minValue = ITL_util<SCALAR>::Min( (SCALAR*)scalarField->getDataFull(), scalarField->getSize() );
		maxValue = ITL_util<SCALAR>::Max( (SCALAR*)scalarField->getDataFull(), scalarField->getSize() );
		// ******************* Need to resolve typecast issue end *************************
	}
	else
	{
		minValue = histMin;
		maxValue = histMax;
	}

	histMin = minValue;
	histMax = maxValue;

	printf("max:%g\tmin:%g",maxValue,minValue);
	h = ( histMax - histMin ) / (float) nBinArray[nIter-1];
	printf("h:%g\n",h);

	computeHistogramBinField_h( scalarField, fieldType, nBinArray[nIter-1], NULL );

	freqListPtr = new double[nBinArray[nIter-1]];

	for( int j=0; j<nBinArray[nIter-1]; j++ )
		freqListPtr[j] = 0;

	//computeFrequencies_h( this->binDatas->grid->nVertices, this->binDatas->datastore->array, freqListPtr );
	computeFrequencies_h( this->binDatas->getSize(), this->binDatas->getDataFull(), freqListPtr );

	freqListPtr2 = new double*[nIter];
	for(int k=0;k<nIter;k++)
		freqListPtr2[k] = new double[nBinArray[nIter-1]];


	for( int j=0; j<nBinArray[nIter-1]; j++ )
	{
		freqListPtr[j] /= (float)N;
		freqListPtr2[nIter-1][j] = freqListPtr[j];
		assert( freqListPtr[j] >= 0 );
		sumOfSquares += ( freqListPtr[j] * freqListPtr[j] );

	}

	firstPart =  2.0f / ((N-1)*h);
	secondPart = (double)( N+1 ) / (double)(( N-1 )*h);
	//secondPart = (double)( N+1 ) / (double)(( N-1 )*N*N*h);
	scoreList[nIter-1] = firstPart - ( secondPart * sumOfSquares );
	//fprintf(crossV,"%d\t %g\t %g\n", nBinArray[nIter-1], h, scoreList[nIter-1]);
	printf( "Number of bins:%d Binwidth:%g Score:%g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );
    //fprintf( crossV, "%d\t %lf\t %lf\n", nBinArray[nIter-1], h, scoreList[nIter-1] );

	//fprintf(crossV,"%d\t1st:%g\t 2nd:%g\t sos:%g\t score:%g\n", nBinArray[nIter-1],firstPart, secondPart, sumOfSquares, scoreList[nIter-1]);
	if( scoreList[nIter-1] < minScore )
	{
		minScore = scoreList[nIter-1];
		minScoreIndex = nIter-1;
	}

	//printf( "==========%d, %g, %g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );

	for( int i = nIter-2; i>=0; i-- )
	{

		h = ( maxValue - minValue ) / (double) nBinArray[i];
		sumOfSquares = 0;
		for( int j=1; j<=nBinArray[i]; j++ )
		{
			freqListPtr2[i][j-1] = freqListPtr2[i+1][j*2-1] + freqListPtr2[i+1][j*2-2];
			assert( freqListPtr2[i][j-1] >= 0 );
			sumOfSquares += ( freqListPtr2[i][j-1] * freqListPtr2[i][j-1] );
		}


		firstPart =  2.0f / ((N-1)*h);
		secondPart = (double)( N+1 ) / (double)(( N-1 )*h);
		//secondPart = (double)( N+1 ) / (double)(( N-1 )*N*N*h);
		scoreList[i] = firstPart - ( secondPart * sumOfSquares );

		//printf( "Number of bins:%d Binwidth:%g Score:%g\n", nBinArray[i], h, scoreList[i] );
		//fprintf( crossV, "%d\t %lf\t %lf\n", nBinArray[i], h, scoreList[i] );
		//fprintf(crossV,"%d\t1st:%g\t 2nd:%g\t sos:%g\t score:%g\n", nBinArray[i],firstPart, secondPart, sumOfSquares, scoreList[i]);

		if( scoreList[i] < minScore )
		{
			minScore = scoreList[i];
			minScoreIndex = i;
		}


	}

/*
	int *freq;
	freq = new int[nBinArray[minScoreIndex]];
	for(int i=0;i<nBinArray[minScoreIndex];i++)
	{
		freq[i] = freqListPtr2[minScoreIndex][i];
	}*/
	//printf("new_start:%d\tnew_Iter:%d",new_start,new_Iter);

	printf( "Cross validation result: optimal number of bins: %d\n", nBinArray[minScoreIndex] );
	cross_nBin=nBinArray[minScoreIndex];
	//printf("minScoreIndex:%d\t nBinArray[minScoreIndex]:%d",minScoreIndex,nBinArray[minScoreIndex]);

/*
	int new_start,new_Iter;
	if(minScoreIndex!=(nIter-1))
	{
		new_start=nBinArray[minScoreIndex];
	    new_Iter = ceil((nBinArray[minScoreIndex]*2 - new_start)/(float)step);
	    crossValidate( "scalar", new_Iter, new_start, step );
	}

	if(minScoreIndex!=0)
	{
		new_start=nBinArray[minScoreIndex]/2;
	    new_Iter = ceil((nBinArray[minScoreIndex]- new_start)/(float)step);
	    crossValidate( "scalar", new_Iter, new_start, step );
	}
*/

	//fclose( crossV );

	delete [] freqListPtr;
	for(int kk=0;kk<nIter;kk++)
		delete [] freqListPtr2[kk];
	delete [] freqListPtr2;
	delete [] scoreList;
	delete [] nBinArray;

	return cross_nBin;

	//for( int i=0; i<nIter; i++ )
	//	delete [] freqListPtr[i];
	//delete [] freqListPtr;

}// End function


void
ITL_histogram::computeHistogramBinField_h( ITL_field_regular<SCALAR> *scalarField, char *fieldType, int nBin, char* binMapFile)
{
	if( strcmp( fieldType, "scalar" ) == 0 )
	{
		// Determine number of bins, if not specified already
		//if( nBin == 0 )		nBin = (int) floor( scalarField->grid->nVertices / 10.0f );
		if( nBin == 0 )		nBin = (int) floor( scalarField->getSize() / 10.0f );

		computeHistogramBinField_Scalar_h( scalarField, nBin );
	}
	/*else if( strcmp( fieldType, "vector" ) == 0 )
	{
		// Determine number of bins, if not specified already
		if( nBin == 0 )		nBin = 360;

		computeHistogramBinField_Vector( nBin, 0, binMapFile );
	}
	else if( strcmp(fieldType, "vector2" ) == 0 )
	{
		// Determine number of bins, if not specified already
		if( nBin == 0 )		nBin = 90;

		computeHistogramBinField_Vector2( nBin );

	}*/
}// End function


void
ITL_histogram::computeHistogramBinField_Scalar_h( ITL_field_regular<SCALAR> *scalarField, int nBin )
{
	//assert( scalarField->datastore->array != NULL );
	assert( scalarField->getDataFull() != NULL );
        SCALAR nextV, minValue, maxValue, rangeValue;

	// The histogram field is padded, pad length is same as neighborhood size of vector field
        //int* lPadHisto = new int[this->grid->nDim];
        //int* hPadHisto = new int[this->grid->nDim];
        //ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
        //ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

	// Initialize the padded scalar field for histogram bins
	if( this->binDatas == NULL )
		//this->binDatas = new ITL_field_regular<int>( scalarField->grid->nDim, scalarField->grid->low, scalarField->grid->high, scalarField->grid->lowPad, scalarField->grid->highPad, scalarField->grid->neighborhoodSizeArray );
	{
		float low[4];
		float high[4];
		scalarField->getBounds( low, high );
		this->binDatas = new ITL_field_regular<int>( scalarField->getNumDim(), low, high );
	}

	// Compute the range over which histogram computation needs to be done
	if( hist_RangeSet == false )
	{
		// Get min-max values of the scalar field
		// ******************* Need to resolve typecast issue begin *************************
		//minValue = ITL_util<SCALAR>::Min( (SCALAR*)scalarField->datastore->array, scalarField->grid->nVertices );
		//maxValue = ITL_util<SCALAR>::Max( (SCALAR*)scalarField->datastore->array, scalarField->grid->nVertices );
		minValue = ITL_util<SCALAR>::Min( (SCALAR*)scalarField->getDataFull(), scalarField->getSize() );
		maxValue = ITL_util<SCALAR>::Max( (SCALAR*)scalarField->getDataFull(), scalarField->getSize() );
		// ******************* Need to resolve typecast issue end *************************
	}
	else
	{
		minValue = histMin;
		maxValue = histMax;
	}
	rangeValue = maxValue - minValue;

	// Compute bin width
	float binWidth = rangeValue / (float)nBin;

	#ifdef DEBUG_MODE
	printf( "Min: %g Max: %g Range: %g of the scalar values\n", minValue, maxValue, rangeValue );
	printf( "Binwidth: %g\n", binWidth );
	#endif

	// Scan through each point of the histogram field
	// and convert field value to bin ID
	int index1d = 0;
	int binId = 0;
	int dimWithPad[4];
	scalarField->getSizeWithPad( dimWithPad );
	for( int z=0; z<dimWithPad[2]; z++ )
	{
		for( int y=0; y<dimWithPad[1]; y++ )
		{
			for( int x=0; x<dimWithPad[0]; x++ )
			{
				// Get scalar value at location
				// ******************* Need to resolve typecast issue *************************
				//nextV = scalarField->datastore->array[index1d];
				nextV = scalarField->getDataAt( index1d );

				// Obtain the binID corresponding to the value at this location
				binId = (int)floor( ( nextV - minValue ) / binWidth  );
				binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
				this->binDatas->setDataAt( index1d, binId );

				// increment to the next grid vertex
				index1d += 1;
			}
		}
	}

    // delete lPadHisto;
    // delete hPadHisto;

}// end function


void
ITL_histogram::setHistogramRange_h( SCALAR minR, SCALAR maxR )
{
	histMin = minR;
	histMax = maxR;
	hist_RangeSet = true;
}


void
ITL_histogram::computeFrequencies_h( int nPoint, int *binIds, double *freqArray )
{
	// Scan through bin Ids and keep count
	for( int i=0; i<nPoint; i++ )
		freqArray[ binIds[i] ] ++;
}

