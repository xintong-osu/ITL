/*
 * ITL_histogrammapper.h
 *
 *  Created on: Apr 8, 2012
 *      Author: abon
 */

#ifndef ITL_HISTOGRAMMAPPER_H_
#define ITL_HISTOGRAMMAPPER_H_

#include "ITL_util.h"
#include "ITL_histogram.h"
#include "ITL_field_regular.h"

template <class T>
class ITL_histogrammapper
{
	ITL_histogram* histogram;
	bool histogramRangeSet;
	T histogramMin, histogramMax;
	T histogramMinArray[2], histogramMaxArray[2];

public:

	ITL_histogrammapper( ITL_histogram* hist )
	{
		this->histogram = hist;
		histogramRangeSet = false;
	}

	void
	setHistogramRange( T m, T M )
	{
		histogramMin = m;
		histogramMax = M;
		histogramRangeSet = true;
	}

	/**
	 * Joint histogram range set function.
	 * Sets the range for a 2D joint histogram.
	 * @param m1 Lower limit of the range of the fist histogram
	 * @param M1 Upper limit of the range of the fist histogram
	 * @param m2 Lower limit of the range of the second histogram
	 * @param m2 Upper limit of the range of the second histogram
	 */
	void
	setJointHistogramRange( T m1, T M1, T m2, T M2 )
	{
		histogramMinArray[0] = m1;
		histogramMaxArray[0] = M1;
		histogramMinArray[1] = m2;
		histogramMaxArray[1] = M2;
		histogramRangeSet = true;
	}

	/**
	 * Histogram bin assignment function for scalar fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Scalar( ITL_field_regular<T>* dataField,
									 ITL_field_regular<int>** binField,
									 int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];

		//assert( this->dataField->datastore->array != NULL );
		assert( dataField->getDataFull() != NULL );
	    SCALAR nextV, rangeValue;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );
			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of the scalar field
			// ******************* Need to resolve typecast issue begin *************************
			//minValue = ITL_util<SCALAR>::Min( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			//maxValue = ITL_util<SCALAR>::Max( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			histogramMin = ITL_util<SCALAR>::Min( (SCALAR*)dataField->getDataFull(), dataField->getSize() );
			histogramMax = ITL_util<SCALAR>::Max( (SCALAR*)dataField->getDataFull(), dataField->getSize() );

			// ******************* Need to resolve typecast issue end *************************
		}
		rangeValue = histogramMax - histogramMin;

		// Compute bin width
		float binWidth = rangeValue / (float)nBin;

		//#ifdef DEBUG_MODE
		printf( "Min: %g Max: %g Range: %g of the scalar values\n", histogramMin, histogramMax, rangeValue );
		printf( "Binwidth: %g\n", binWidth );
		//#endif

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		int dimWithPad[4];
		(*binField)->getSizeWithPad( dimWithPad );
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get scalar value at location
					// ******************* Need to resolve typecast issue *************************
					nextV = dataField->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV - histogramMin ) / binWidth  );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binField)->setDataAt( index1d, binId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

	}// end function


	/**
	 * Histogram bin assignment function for vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Vector( ITL_field_regular<T>* dataField,
			 	 	 	 	 	 	 ITL_field_regular<int>** binField,
			 	 	 	 	 	 	 int nBin, int iRes = 0, char* binMapFile = NULL )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField->getDataFull() != NULL );
		VECTOR3 nextV;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		(*binField)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", this->dataField->grid->dimWithPad[0], this->dataField->grid->dimWithPad[1], this->dataField->grid->dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					nextV = (VECTOR3)dataField->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location
					(*binField)->setDataAt( index1d, histogram->get_bin_number_3D( nextV, iRes ) );

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
	void
	computeHistogramBinField_Vector2( ITL_field_regular<T>* dataField,
	   	 	 	 	 	 	 	      ITL_field_regular<int>** binField,
	 	 	 	 	 	 	 	 	  int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField->getDataFull() != NULL );
		T *nextV = new T();

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;

		(*binField)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
				//printf( "%d %d %d %d %d %d %d\n", this->dataField->grid->dimWithPad[0], this->dataField->grid->dimWithPad[1], this->dataField->grid->dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					*nextV = dataField->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location
					(*binField)->setDataAt( index1d, histogram->get_bin_number_2D( *nextV, nBin ) );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

        delete nextV;

	}// end function

	void
	computeJointHistogramBinField_Scalar( ITL_field_regular<T>* dataField1,
	  								      ITL_field_regular<T>* dataField2,
	 	 	 	   	   	   	   	   	   	  ITL_field_regular<int>** binField,
	 	 	 	   	   	   	   	   	   	  int nBin )
	{
		cout << "Scalar" << endl;
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		T nextVField1, nextVField2;

		//assert( this->dataField1->datastore->array != NULL );
		//assert( this->dataField2->datastore->array != NULL );
		assert( dataField1->getDataFull() != NULL );
		assert( dataField2->getDataFull() != NULL );

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField1->getBounds( low, high );
			dataField1->getPadSize( lowPad, highPad );
			dataField1->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField1->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of both the scalar fields
			histogramMinArray[0] = ITL_util<T>::Min( dataField1->getDataFull(), dataField1->getSize() );
			histogramMaxArray[0] = ITL_util<T>::Max( dataField1->getDataFull(), dataField1->getSize() );
			histogramMinArray[1] = ITL_util<T>::Min( dataField2->getDataFull(), dataField2->getSize() );
			histogramMaxArray[1] = ITL_util<T>::Max( dataField2->getDataFull(), dataField2->getSize() );
		}

		// Compute bin widths along each dimension
		T rangeValue1 = histogramMaxArray[0] - histogramMinArray[0];
		float binWidth1 = rangeValue1 / (float)nBin;
		T rangeValue2 = histogramMaxArray[1] - histogramMinArray[1];
		float binWidth2 = rangeValue2 / (float)nBin;

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId1, binId2, combinedBinId;
		int dimWithPad[4];
		dataField1->getSizeWithPad( dimWithPad );
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					nextVField1 = dataField1->getDataAt( index1d );
					nextVField2 = dataField1->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location from both fields individually
					binId1 = (int)floor( ( nextVField1 - histogramMinArray[0] ) / binWidth1  );
					binId1 = ITL_util<int>::clamp( binId1, 0, nBin-1 );
					binId2 = (int)floor( ( nextVField2 - histogramMinArray[1] ) / binWidth2  );
					binId2 = ITL_util<int>::clamp( binId2, 0, nBin-1 );

					// Combine the two bin indices from the two fields
					combinedBinId = binId2 * nBin + binId1;

					// Store the bin ID
					(*binField)->setDataAt( index1d, combinedBinId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

		// delete lPadHisto;
		// delete hPadHisto;

	}// end function

	/**
	 * Joint histogram bin assignment function for vector field.
	 * @param nBin Number of bins to be used in joint histogram computation.
	 */
	void
	computeJointHistogramBinField_Vector( ITL_field_regular<T>* dataField1,
										  ITL_field_regular<T>* dataField2,
										  ITL_field_regular<int>** binField,
										  int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField1->getDataFull() != NULL );
		assert( dataField2->getDataFull() != NULL );
		VECTOR3 nextVField1;
		VECTOR3 nextVField2;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField1->getBounds( low, high );
			dataField1->getPadSize( lowPad, highPad );
			dataField1->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField1->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		dataField1->getSizeWithPad( dimWithPad );
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					nextVField1 = dataField1->getDataAt( index1d );
					nextVField2 = dataField2->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location from both fields individually
					int binValue1 = histogram->get_bin_number_3D( nextVField1 );
					int binValue2 = histogram->get_bin_number_3D( nextVField2 );

					// Combine the two bin indices from the two fields
					int combinedBin = binValue2 * nBin + binValue1;

					// Store the bin ID
					(*binField)->setDataAt( index1d, combinedBin );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

	}// end function

	static void
	computeHistogramFrequencies( ITL_field_regular<int>** binField,
								 int* freqList,
								 int nBin )
	{
		assert( binField != NULL );
		assert( freqList != NULL );

		// Scan through bin Ids and keep count
		for( int i=0; i<nBin; i++ )
			freqList[i] = 0;

		printf( "here\n");
		for( int i=0; i<(*binField)->getSize(); i++ )
		{
			//cout << (*binField)->getDataAt(i) << endl;
			freqList[ (*binField)->getDataAt(i) ] ++;
		}

	}

	static void
	computeHistogramFrequencies( int *binIds, float* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
	}

	static void
	computeHistogramFrequencies( int *binIds, double* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
	}

	/**
	  * Cross-validation function.
	  */
	void crossValidate( ITL_field_regular<T>* dataField,
						char *fieldType,
						int nIter,
						int start, int step )
	{
		ITL_field_regular<int>* binField = NULL;


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
		int N = dataField->getSize();
		cout << "N " << N << endl;
		double firstPart = 0;
		double secondPart = 0;
		double sumOfSquares = 0;


		//FILE *histofile = fopen( "histo.csv", "w" );
		for( int i = 0; i<nIter; i++ )
		{
			h = ( histogramMax - histogramMin ) / (float) nBinArray[i];

			if( binField != NULL )
				delete binField;
			if( strcmp( fieldType, "scalar" ) == 0 )
				computeHistogramBinField_Scalar( dataField, &binField, nBinArray[nIter-1] );
			else if( strcmp( fieldType, "vector" ) == 0 )
				computeHistogramBinField_Vector( dataField, &binField, nBinArray[nIter-1] );
			else if( strcmp( fieldType, "vector2" ) == 0 )
				computeHistogramBinField_Vector2( dataField, &binField, nBinArray[nIter-1] );

			freqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				freqListPtr[j] = 0;

			//computeFrequencies( this->binData->grid->nVertices, this->binData->datastore->array, freqListPtr );
			computeHistogramFrequencies( binField->getDataFull(),
					 	 	 	 	 	 freqListPtr,
										 binField->getSize() );

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

		//fclose( histofile );

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
	int crossValidateSpeedUp( ITL_field_regular<T>* dataField,
							  char *fieldType,
							  int nMax,
							  int start, int step )
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
		int N = dataField->getSize();

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
			histogramMin = ITL_util<SCALAR>::Min( (SCALAR*)dataField->getDataFull(), dataField->getSize() );
			histogramMax = ITL_util<SCALAR>::Max( (SCALAR*)dataField->getDataFull(), dataField->getSize() );
			// ******************* Need to resolve typecast issue end *************************
		}

		printf("max:%g\tmin:%g",maxValue,minValue);
		h = ( histogramMax - histogramMin ) / (float) nBinArray[nIter-1];
		printf("h:%g\n",h);

		ITL_field_regular<int>* binField = NULL;
		if( strcmp( fieldType, "scalar" ) == 0 )
			computeHistogramBinField_Scalar( dataField, &binField, nBinArray[nIter-1] );
		else if( strcmp( fieldType, "vector" ) == 0 )
			computeHistogramBinField_Vector( dataField, &binField, nBinArray[nIter-1] );
		else if( strcmp( fieldType, "vector2" ) == 0 )
			computeHistogramBinField_Vector2( dataField, &binField, nBinArray[nIter-1] );

		freqListPtr = new double[nBinArray[nIter-1]];

		for( int j=0; j<nBinArray[nIter-1]; j++ )
			freqListPtr[j] = 0;

		//computeFrequencies( this->binData->grid->nVertices, this->binData->datastore->array, freqListPtr );
		ITL_histogrammapper::computeHistogramFrequencies( binField->getDataFull(), freqListPtr,  binField->getSize() );

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

};

#endif
/* ITL_HISTOGRAMMAPPER_H_ */
