/*
 * ITL_histogrammapper.h
 *
 *  Created on: Apr 8, 2012
 *  Author: abon
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
	SCALAR histogramMin, histogramMax;
	SCALAR histogramMinArray[2], histogramMaxArray[2];

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
		SCALAR nextV, rangeValue;

		//cout << "in" << endl;
		assert( dataField->getDataFull() != NULL );
		//cout << "in2" << endl;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}
		//cout << "in3" << endl;

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			//cout << "setting histogram range" << endl;
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

		#ifdef DEBUG_MODE
		printf( "Min: %g Max: %g Range: %g of the scalar values\n", histogramMin, histogramMax, rangeValue );
		printf( "Binwidth: %g\n", binWidth );
		#endif

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		int dimWithPad[4];
		(*binField)->getSizeWithPad( dimWithPad );
		//cout << "starting iteration" << endl;
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
		///cout << "data scan completed" << endl;

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

		for( int i=0; i<(*binField)->getSize(); i++ )
		{
			//cout << (*binField)->getDataAt(i) << endl;
			freqList[ (*binField)->getDataAt(i) ] ++;
		}

	}

	static void
	computeHistogramFrequencies( int *binIds, int* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
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
						int start, int step,
						int *nBinArray,
						double *scoreList,
						int* nBinOptimal )
	{
		int nDim = dataField->getNumDim();
		float lowF[3], highF[3];
		dataField->getBounds( lowF, highF );
		int N = (double)dataField->getSize();
		int *freqListPtr = NULL;
		double *normFreqListPtr = NULL;

		// Create bin field
		ITL_field_regular<int>* binField = new ITL_field_regular<int>( nDim, lowF, highF );

		// Compute number of bins for different iterations
		for( int i = 0; i<nIter; i++ )
		{
			nBinArray[i] = start + i*step;
		}

		// Iterate
		double minScore = 10000000;
		int minScoreIndex = -1;
		double h = 0;

		// Scan through different number of bins
		for( int i = 0; i<nIter; i++ )
		{
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			// Compute bin field
			if( strcmp( fieldType, "scalar" ) == 0 )
				computeHistogramBinField_Scalar( dataField, &binField, nBinArray[i] );
			else if( strcmp( fieldType, "vector" ) == 0 )
				computeHistogramBinField_Vector( dataField, &binField, nBinArray[i] );
			else if( strcmp( fieldType, "vector2" ) == 0 )
				computeHistogramBinField_Vector2( dataField, &binField, nBinArray[i] );

			// Compute frequencies
			freqListPtr = new int[nBinArray[i]];
			memset( freqListPtr, 0, sizeof(int)*nBinArray[i] );
			computeHistogramFrequencies( binField->getDataFull(),
					 	 	 	 	 	 freqListPtr,
										 binField->getSize() );

			// Normalize frequencies
			normFreqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr[j] / (double)N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			//#endif

			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			delete [] freqListPtr;
			delete [] normFreqListPtr;

		}// end for : iterations

		delete binField;

		(*nBinOptimal) = nBinArray[minScoreIndex];

	}// End function


	/** added by Tzu-Hsuan
	 *
	 */
	void
	crossValidate_optimized( ITL_field_regular<T>* dataField,
	 					  	  char *fieldType,
							  int nIter, int nBinMax,
  							  int *nBinArray,
							  double* scoreList,
							  int* nBinOptimal )
	{
		int nDim = dataField->getNumDim();
		float lowF[3], highF[3];
		dataField->getBounds( lowF, highF );
		int N = dataField->getSize();
		int* freqListPtr = NULL;
		int* freqListPtr2 = NULL;
		double* normFreqListPtr = NULL;
		int cross_nBin;
		double h = 0;

		// Create bin field
		ITL_field_regular<int>* binField = new ITL_field_regular<int>( nDim, lowF, highF );

		// Compute number of bins for different iterations
		for( int i = 0; i<nIter; i++ )
		{
			nBinArray[nIter-i-1] = (int)(nBinMax/(pow(2.0,(float)i)));

		}

		// Compute bin field from data only in the first iteration
		if( strcmp( fieldType, "scalar" ) == 0 )
			computeHistogramBinField_Scalar( dataField, &binField, nBinArray[nIter-1] );
		else if( strcmp( fieldType, "vector" ) == 0 )
			computeHistogramBinField_Vector( dataField, &binField, nBinArray[nIter-1] );
		else if( strcmp( fieldType, "vector2" ) == 0 )
			computeHistogramBinField_Vector2( dataField, &binField, nBinArray[nIter-1] );

		// Compute frequencies
		freqListPtr = new int[nBinArray[nIter-1]];
		memset( freqListPtr, 0, sizeof(int)*nBinArray[nIter-1] );
		ITL_histogrammapper::computeHistogramFrequencies( binField->getDataFull(), freqListPtr,  binField->getSize() );
		h = ( histogramMax - histogramMin ) / (double) nBinArray[nIter-1];

		// Normalize frequencies
		normFreqListPtr = new double[nBinArray[nIter-1]];
		for( int j=0; j<nBinArray[nIter-1]; j++ )
			normFreqListPtr[j] =  freqListPtr[j] / (double)N;

		// Compute cross-validation score
		scoreList[nIter-1] = computeCrossValidationScore( normFreqListPtr, nBinArray[nIter-1], N, h );

		#ifdef DEBUG_MODE
		printf( "%d, %g, %g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );
		#endif

		// Record latest minimum score
		double minScore = scoreList[nIter-1];
		int minScoreIndex = (nIter-1);

		// Clear temp resources
		delete [] normFreqListPtr;

		// Update
		for( int i = nIter-2; i>=0; i-- )
		{
			// Update bin width
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			// Update frequencies
			freqListPtr2 = new int[nBinArray[i]];
			for( int j=0; j<=nBinArray[i]; j++ )
				freqListPtr2[j] = freqListPtr[j*2] + freqListPtr[j*2+1];

			// Normalize updated frequencies
			normFreqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr2[j] / (double)N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			//#endif

			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			// Store the current frequent list for next iteration
			memcpy( freqListPtr, freqListPtr2, sizeof(int)*nBinArray[i] );

			delete [] normFreqListPtr;
			delete [] freqListPtr2;

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];

		delete [] freqListPtr;

	}// End function

	void
	crossValidate_optimized2( int* freqList,
							  int nIter, int nBinMax, int N,
							  int *nBinArray,
  							  double* scoreList,
  							  int* nBinOptimal )
	{
		int* freqListPtr = new int[nBinMax];
		memcpy( freqListPtr, freqList, sizeof(int)*nBinMax );
		int* freqListPtr2 = NULL;
		double* normFreqListPtr = NULL;
		double h = 0;

		// Compute number of bins for different iterations
		for( int i = 0; i<nIter; i++ )
			nBinArray[nIter-i-1] = (int)(nBinMax/(pow(2.0,(float)i)));

		// Normalize frequencies
		normFreqListPtr = new double[nBinArray[nIter-1]];
		for( int j=0; j<nBinArray[nIter-1]; j++ )
			normFreqListPtr[j] =  freqListPtr[j] / (double)N;

		// Compute cross-validation score
		h = ( histogramMax - histogramMin ) / (double) nBinArray[nIter-1];
		scoreList[nIter-1] = computeCrossValidationScore( normFreqListPtr, nBinArray[nIter-1], N, h );

		#ifdef DEBUG_MODE
		printf( "%d, %g, %g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );
		#endif

		// Record latest minimum score
		double minScore = scoreList[nIter-1];
		int minScoreIndex = (nIter-1);

		// Clear temp resources
		delete [] normFreqListPtr;

		// Update
		for( int i = nIter-2; i>=0; i-- )
		{
			// Update bin width
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			// Update frequencies
			freqListPtr2 = new int[nBinArray[i]];
			for( int j=0; j<=nBinArray[i]; j++ )
				freqListPtr2[j] = freqListPtr[j*2] + freqListPtr[j*2+1];

			// Normalize updated frequencies
			normFreqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr2[j] / (double)N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			#endif

			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			// Store the current frequent list for next iteration
			memcpy( freqListPtr, freqListPtr2, sizeof(int)*nBinArray[i] );

			delete [] normFreqListPtr;
			delete [] freqListPtr2;

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];

		delete [] freqListPtr;

	}// End function

	/**
	 * Cross-validation score computation function.
	 *
	 */
	double
	computeCrossValidationScore( double* fList, int nBin, int N, double h )
	{
		double sumOfSquares = 0;
		for( int j = 0; j<nBin; j++ )
		{
			assert( fList[j] >= 0 );
			sumOfSquares += ( fList[j] * fList[j] );
		}

		double firstPart =  2.0f / ((N-1)*h);
		//double secondPart = (double)( N+1 ) / (double)(( N-1 )*h);
		double secondPart = (double)( N+1 ) / (double)(( N-1 )*N*N*h);
		double score = firstPart - ( secondPart * sumOfSquares );

		return score;

	}// end function

};

#endif
/* ITL_HISTOGRAMMAPPER_H_ */
