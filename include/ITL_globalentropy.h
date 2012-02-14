/**
 * Gloabl entropy computation class.
 * Created on: April 02, 2011.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GLOBALENTROPY_H_
#define ITL_GLOBALENTROPY_H_

#include "ITL_header.h"
#include "ITL_field_regular.h"
#include "ITL_histogram.h"
#include "ITL_entropycore.h"

template <class T>
class ITL_globalentropy
{
public:

	ITL_field_regular<T> *dataField;		/**< A scalar field containing field data */
	ITL_field_regular<int>* binData;		/**< A scalar field containing histogram bins corresponding to field points. */
	int *freqList;

	T histogramMin;
	T histogramMax;
	bool histogramRangeSet;

	float globalEntropy;				/**< Value of computed global entropy of the field. */
	float* probarray;   				/**< Probability array used in entropy computation. */
	
	ITL_histogram *histogram;			// ADD-BY-ABON 11/07/2011

public:

	/**
	 * Constructor Type 1.
	 */
	ITL_globalentropy( ITL_field_regular<T> *f, ITL_histogram *hist )
	{
		dataField = f;
		binData = NULL;
		freqList = NULL;
		histogramRangeSet = false;
		histogram = hist;

	}// End constructor

	/**
	 * Constructor Type 2.
	 */
	ITL_globalentropy( ITL_field_regular<T> *dataF, ITL_field_regular<int> *binF, ITL_histogram *hist  )
	{
		dataField = dataF;
		binData = binF;
		freqList = NULL;
		histogramRangeSet = false;
		histogram = hist;

	}// End constructor

	/**
	 * Histogram bin assignment function.
	 * Calls the appropriate function based on field type
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField( char *fieldType, int nBin = 0, char* binMapFile = NULL )
	{
		if( strcmp( fieldType, "scalar" ) == 0 )
		{
			// Determine number of bins, if not specified already
			if( nBin == 0 )		nBin = (int) floor( this->dataField->grid->nVertices / 10.0f );
			
			computeHistogramBinField_Scalar( nBin );
		}
		else if( strcmp( fieldType, "vector" ) == 0 )
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

		}
	}// End function
	
	/**
	 * Histogram bin assignment function for scalar fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField_Scalar( int nBin )
	{
		assert( this->dataField->datastore->array != NULL );
	    SCALAR nextV, minValue, maxValue, rangeValue;
		
		// The histogram field is padded, pad length is same as neighborhood size of vector field
	        //int* lPadHisto = new int[this->grid->nDim];
	        //int* hPadHisto = new int[this->grid->nDim];
	        //ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
	        //ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

		// Initialize the padded scalar field for histogram bins
		if( this->binData == NULL )
			this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim, this->dataField->grid->low, this->dataField->grid->high, this->dataField->grid->lowPad, this->dataField->grid->highPad, this->dataField->grid->neighborhoodSizeArray );
					
		
		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of the scalar field
			// ******************* Need to resolve typecast issue begin *************************
			minValue = ITL_util<SCALAR>::Min( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			maxValue = ITL_util<SCALAR>::Max( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			// ******************* Need to resolve typecast issue end *************************
		}
		else
		{
			minValue = histogramMin;
			maxValue = histogramMax;
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
		for( int z=0; z<this->dataField->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField->grid->dimWithPad[0]; x++ )
				{
					// Get scalar value at location
					// ******************* Need to resolve typecast issue *************************
					nextV = this->dataField->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV - minValue ) / binWidth  );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					this->binData->setDataAt( index1d, binId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

        // delete lPadHisto;
        // delete hPadHisto;

	}// end function


	/**
	 * Histogram bin assignment function for vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField_Vector( int nBin, int iRes = 0, char* binMapFile = NULL )
	{
		assert( this->dataField->datastore->array != NULL );
		VECTOR3 nextV;

		// The histogram field is padded, pad length is same as neighborhood size of vector field
		//int* lPadHisto = new int[this->grid->nDim];
		//int* hPadHisto = new int[this->grid->nDim];
		//ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
		//ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

		// Initialize the padded scalar field for histogram bins
		if( this->binData == NULL )
			//this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim,
			//											this->dataField->grid->low, this->dataField->grid->high,
			//											this->dataField->grid->lowPad, this->dataField->grid->highPad,
			//											//lPadHisto, hPadHisto,
			//											this->dataField->grid->neighborhoodSize );
			this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim,
														this->dataField->grid->low, this->dataField->grid->high,
														this->dataField->grid->lowPad, this->dataField->grid->highPad,
														//lPadHisto, hPadHisto,
														this->dataField->grid->neighborhoodSizeArray );

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		//printf( "nVerts: %d\n", this->dataField->grid->nVertices );
		for( int z=0; z<this->dataField->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField->grid->dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", this->dataField->grid->dimWithPad[0], this->dataField->grid->dimWithPad[1], this->dataField->grid->dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					nextV = (VECTOR3)this->dataField->datastore->array[index1d];
										
					// Obtain the binID corresponding to the value at this location
					this->binData->setDataAt( index1d, histogram->get_bin_number_3D( nextV, iRes ) );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

        // delete lPadHisto;
        // delete hPadHisto;

	}// end function

	/**
	 * Histogram bin assignment function for 2D vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField_Vector2( int nBin )
	{
		assert( this->dataField->datastore->array != NULL );
		T *nextV = new T();

		// The histogram field is padded, pad length is same as neighborhood size of vector field
		//int* lPadHisto = new int[this->grid->nDim];
		//int* hPadHisto = new int[this->grid->nDim];
		//ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
		//ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

		// Initialize the padded scalar field for histogram bins
		if( this->binData == NULL )
			//this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim,
			//											this->dataField->grid->low, this->dataField->grid->high,
			//											this->dataField->grid->lowPad, this->dataField->grid->highPad,
			//											//lPadHisto, hPadHisto,
			//											this->dataField->grid->neighborhoodSize );
			this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim,
														this->dataField->grid->low, this->dataField->grid->high,
														this->dataField->grid->lowPad, this->dataField->grid->highPad,
														//lPadHisto, hPadHisto,
														this->dataField->grid->neighborhoodSizeArray );

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		//printf( "nVerts: %d\n", this->dataField->grid->nVertices );
		for( int z=0; z<this->dataField->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField->grid->dimWithPad[0]; x++ )
				{
				//printf( "%d %d %d %d %d %d %d\n", this->dataField->grid->dimWithPad[0], this->dataField->grid->dimWithPad[1], this->dataField->grid->dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					*nextV = this->dataField->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location
					this->binData->setDataAt( index1d, histogram->get_bin_number_2D( *nextV, nBin ) );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

        delete nextV;
        // delete lPadHisto;
        // delete hPadHisto;

	}// end function


	/**
	 * Entropy computation function.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBins Number of bins used in histogram computation/Number of sample points in KDE estimation.
	 * @paran toNormalize boolean flag, TRUE indicates that the computed entropy value needs to be normalized
	 */
	void computeGlobalEntropyOfField( int nBin, bool toNormalize, int method = 0 )
	{
		if( method == 0 )
		{
			// Check if frequencies are already computed from bin Ids
			// If not, only then compute frequencies from bin Ids
			if( this->freqList == NULL )
				freqList = new int[nBin];
			computeHistogramFrequencies( nBin );
			
			// Compute entropy from frequency list			
			globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( freqList,
										binData->grid->nVertices, nBin, toNormalize );	
	
		}		
		else
			//globalEntropy = ITL_entropycore::computeEntropy_HistogramBased(
			//								this->binData->datastore->array,
			//								this->binData->grid->nVertices,
			//								nBin, toNormalize );		
			globalEntropy = ITL_entropycore::computeEntropy_KDEBased(
											this->dataField->datastore->array,
											this->dataField->grid->nVertices,
											0.0, toNormalize );
	}// End function

	/**
	  * Cross-validation function.
	  */
	void crossValidate( char *fieldType, int nIter, int start, int step )
	{
		// Assume method is 0 (histogram-based)
		// Allocate memory
		double *scoreList = new double[nIter];
		int *nBinArray = new int[nIter];
		//int **freqListPtr = new int*[nIter];
		double *freqListPtr = NULL;
		
		for( int i = 0; i<nIter; i++ )
		{
			nBinArray[i] = start + i*step;
		}

		// Iterate
		double minScore = 10000000;
		int minScoreIndex = -1;
		float h = 0;
		int N = this->dataField->grid->nVertices;
		cout << "N " << N << endl;
		double firstPart = 0;
		double secondPart = 0;
		double sumOfSquares = 0;

		
		FILE *histofile = fopen( "histo.csv", "w" );

		for( int i = 0; i<nIter; i++ )
		{
			h = ( histogramMax - histogramMin ) / (float) nBinArray[i];
						
			computeHistogramBinField( fieldType, nBinArray[i] );

			freqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				freqListPtr[j] = 0;

			computeFrequencies( this->binData->grid->nVertices, this->binData->datastore->array, freqListPtr );

			//for( int j=0; j<nBinArray[i]; j++ )
			//{
			//	if( j == nBinArray[i]-1 ) fprintf( histofile, "%g\n", freqListPtr[j] );
			//	else 			  fprintf( histofile, "%g, ", freqListPtr[j] );
			//}
			

			for( int j=0; j<nBinArray[i]; j++ )
				freqListPtr[j] /= (float)N;

			sumOfSquares = 0;
			for( int j = 0; j<nBinArray[i]; j++ )
			{
				//cout << j << ": " << freqListPtr[j] << endl;
				assert( freqListPtr[j] >= 0 );
				sumOfSquares += ( freqListPtr[j] * freqListPtr[j] );
			}
			//printf( "Histogram range: %g %g\n", histogramMin, histogramMax ); 
				//cout << sumOfSquares << endl;
			
			firstPart =  2.0f / ((N-1)*h);
			secondPart = (double)( N+1 ) / (double)(( N-1 )*h);
			//secondPart = (double)( N+1 ) / (double)(( N-1 )*N*N*h);
			scoreList[i] = firstPart - ( secondPart * sumOfSquares );

			//printf( "Number of bins:%d Binwidth:%g Score:%g\n", nBinArray[i], h, scoreList[i] );
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );


			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			delete [] freqListPtr;
		}

		fclose( histofile );

		printf( "Cross validation result: optimal number of bins: %d\n", nBinArray[minScoreIndex] );	

		delete [] scoreList;
		delete [] nBinArray;
		//for( int i=0; i<nIter; i++ )
		//	delete [] freqListPtr[i];
		//delete [] freqListPtr;

	}// End function

	/** added by Tzu-Hsuan
	 *
	 */
	int crossValidateSpeedUp( char *fieldType, int nMax, int start, int step )
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
			nBinArray[nIter-i-1] = (int)(nMax/(pow(2.0,(float)i)));

		}

		for(int i=0;i<nIter;i++)
			printf("nBinArray[i]:%d\n",nBinArray[i]);

		// Iterate
		double minScore = 10000000;
		int minScoreIndex = -1;
		float h = 0;
		int N = this->dataField->grid->nVertices;

		//FILE *crossV;
		//crossV = fopen("crossV.txt","w+");
		double firstPart = 0;
		double secondPart = 0;
		double sumOfSquares = 0;

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of the scalar field
			// ******************* Need to resolve typecast issue begin *************************
			minValue = ITL_util<SCALAR>::Min( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			maxValue = ITL_util<SCALAR>::Max( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			// ******************* Need to resolve typecast issue end *************************
		}
		else
		{
			minValue = histogramMin;
			maxValue = histogramMax;
		}

		histogramMin = minValue;
		histogramMax = maxValue;

		printf("max:%g\tmin:%g",maxValue,minValue);
		h = ( histogramMax - histogramMin ) / (float) nBinArray[nIter-1];
		printf("h:%g\n",h);

		computeHistogramBinField( fieldType, nBinArray[nIter-1] );

		freqListPtr = new double[nBinArray[nIter-1]];

		for( int j=0; j<nBinArray[nIter-1]; j++ )
			freqListPtr[j] = 0;

		computeFrequencies( this->binData->grid->nVertices, this->binData->datastore->array, freqListPtr );

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
		printf( "Optimal Number of bins:%d  Binwidth:%g  Score:%g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );
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
		}
		*/
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

	void computeHistogramFrequencies( int nBin )
	{
		// Scan through bin Ids and keep count
		if( freqList == NULL ) freqList = new int[nBin];
		for( int i=0; i<nBin; i++ )
			freqList[i] = 0;
		for( int i=0; i<binData->grid->nVertices; i++ )
			freqList[ binData->datastore->array[i] ] ++;
	}

	void computeFrequencies( int nPoint, int *binIds, double *freqArray )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqArray[ binIds[i] ] ++;
	}

	/**
	  * Historgam range set function
	  * Sets the range over which historgam computation is to be done.
	  * @param minR lower limit of the range.
	  * @param maxR upper limit of the range.
	  */
	void setHistogramRange( T minR, T maxR )
	{
		histogramMin = minR;
		histogramMax = maxR;
		histogramRangeSet = true;		
	}

	/**
	 * Histogram frequency data accessor function.
	 * Returns pointer to integer array storing the histogram freqencies.
	 */
	void getHistogramFrequencies( int nBin, int *freqlist )
	{
		assert( this->freqList != NULL );
		assert( freqlist != NULL );

		// Copy frequency for each bin
		memcpy( freqlist, this->freqList, sizeof( int ) * nBin );	

	}// end function

	/**
	 * Histogram frequency data accessor function.
	 * Returns pointer to integer array storing the histogram freqencies.
	 */
	void getHistogramFrequencies( int nBin, float *freqlist )
	{
		assert( this->freqList != NULL );
		assert( freqlist != NULL );
		float *tempFreqList = new float[nBin];
		for( int i=0; i<nBin; i++ )
			tempFreqList[i] = freqList[i] / (float)this->binData->grid->nVertices;

		// Copy frequency for each bin
		memcpy( freqlist, tempFreqList, sizeof( float ) * nBin );

		delete [] tempFreqList;

	}// end function

	/**
	 * Histogram bin field accessor function.
	 * Returns pointer to integer array storing the histogram freqencies.
	 */
	ITL_field_regular<int>* getHistogramBinField()
	{
		return this->binData;
	}// end function

	/**
	 * Entropy field accessor function.
	 * Returns computed entropy field.
	 * @return pointer to entropy field.
	 */
	float getGlobalEntropy()
	{
		return this->globalEntropy;
	}// end function

	/** Default Destructor
	 *
	 */
	~ITL_globalentropy()
	{
		if( this->binData != NULL ) 	delete this->binData;
	}// end destructor
};

#endif
/* ITL_GLOBALENTROPY_H_ */
