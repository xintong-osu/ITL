/**
 * Joint entropy computation class for the entire domain.
 * Created on: Jul 19, 2011.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GLOBALJOINTENTROPY_H_
#define ITL_GLOBALJOINTENTROPY_H_

#include "ITL_header.h"
#include "ITL_field_regular.h"
#include "ITL_entropycore.h"

template <class T>
class ITL_globaljointentropy
{
public:

	ITL_field_regular<T> *dataField1;
	ITL_field_regular<T> *dataField2;
	ITL_field_regular<int> *binData;
	
	int *jointFreqList;
	float globalJointEntropy;

	T histogramMin1, histogramMax1, histogramMin2, histogramMax2;
	bool histogramRangeSet;

	ITL_histogram *histogram;	// ADD-BY-ABON 11/07/2011

public:

	/**
	 * Default Constructor.
	 */
	ITL_globaljointentropy( ITL_field_regular<T> *f1, ITL_field_regular<T> *f2, ITL_histogram *hist )
	{
		this->dataField1 = f1;
		this->dataField2 = f2;
		this->binData = NULL;
		this->jointFreqList = NULL;
		histogramRangeSet = false;
		histogram = hist;
	}
	
	/**
	 * Joint histogram bin assignment function.
	 * Calls the appropriate function based on field type
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeJointHistogramBinField( char *fieldType, int nBin = 0 )
	{
		if( strcmp( fieldType, "scalar" ) == 0 )
		{
			// Determine number of bins, if not specified already
			if( nBin == 0 )		nBin = (int) floor( this->dataField1->grid->nVertices / 10.0f );
			
			computeJointHistogramBinField_Scalar( nBin );
		}
		else if( strcmp( fieldType, "vector" ) == 0 )
		{
			// Determine number of bins, if not specified already
			if( nBin == 0 )		nBin = 360;
			
			computeJointHistogramBinField_Vector( nBin );
		}		
	}// End function
	
	/**
	 * Joint histogram bin assignment function for scalar field.
	 * @param nBin Number of bins to be used in histogram computation.
	 */
	void computeJointHistogramBinField_Scalar( int nBin )
	{
		assert( this->dataField1->datastore->array != NULL );
		assert( this->dataField2->datastore->array != NULL );
		SCALAR nextVField1, nextVField2;
		SCALAR minValue1, maxValue1, minValue2, maxValue2;

		// The histogram field is padded, pad length is same as neighborhood size of vector field
		//int* lPadHisto = new int[this->grid->nDim];
		//int* hPadHisto = new int[this->grid->nDim];
		//ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
		//ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

		// Initialize the padded scalar field for histogram bins
		if( this->binData == NULL )
			this->binData = new ITL_field_regular<int>( this->dataField1->grid->nDim,
													this->dataField1->grid->low, this->dataField1->grid->high,
													this->dataField1->grid->lowPad, this->dataField1->grid->highPad,
													this->dataField1->grid->neighborhoodSizeArray );
		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of both the scalar fields
			minValue1 = ITL_util<T>::Min( this->dataField1->datastore->array, this->dataField1->grid->nVertices );
			maxValue1 = ITL_util<T>::Max( this->dataField1->datastore->array, this->dataField1->grid->nVertices );
			minValue2 = ITL_util<T>::Min( this->dataField2->datastore->array, this->dataField2->grid->nVertices );
			maxValue2 = ITL_util<T>::Max( this->dataField2->datastore->array, this->dataField2->grid->nVertices );
		}
		else
		{
			minValue1 = histogramMin1;
			maxValue1 = histogramMax1;
			minValue2 = histogramMin2;
			maxValue2 = histogramMax2;				
		}
		
		// Compute bin widths along each dimension
		T rangeValue1 = maxValue1 - minValue1;
		float binWidth1 = rangeValue1 / (float)nBin;
		T rangeValue2 = maxValue2 - minValue2;
		float binWidth2 = rangeValue2 / (float)nBin;

		#ifdef DEBUG_MODE
		printf( "Field 1:\n");
		printf( "Min: %g Max: %g Range: %g of the scalar values\n", minValue1, maxValue1, rangeValue1 );
		printf( "Binwidth: %g\n", binWidth1 );
		printf( "Field 2:\n");
		printf( "Min: %g Max: %g Range: %g of the scalar values\n", minValue2, maxValue2, rangeValue2 );
		printf( "Binwidth: %g\n", binWidth2 );
		#endif

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId1, binId2, combinedBinId;
		for( int z=0; z<this->dataField1->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField1->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField1->grid->dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					nextVField1 = this->dataField1->datastore->array[index1d];
					nextVField2 = this->dataField2->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location from both fields individually
					binId1 = (int)floor( ( nextVField1 - minValue1 ) / binWidth1  );
					binId1 = ITL_util<int>::clamp( binId1, 0, nBin-1 );
					binId2 = (int)floor( ( nextVField2 - minValue2 ) / binWidth2  );
					binId2 = ITL_util<int>::clamp( binId2, 0, nBin-1 );
					
					// Combine the two bin indices from the two fields
					combinedBinId = binId2 * nBin + binId1;

					// Store the bin ID
					this->binData->setDataAt( index1d, combinedBinId );

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
	void computeJointHistogramBinField_Vector( int nBin )
	{
		assert( this->dataField1->datastore->array != NULL );
		assert( this->dataField2->datastore->array != NULL );
		T* nextVField1 = new T();
		T* nextVField2 = new T();

		// The histogram field is padded, pad length is same as neighborhood size of vector field
		//int* lPadHisto = new int[this->grid->nDim];
		//int* hPadHisto = new int[this->grid->nDim];
		//ITL_util<int>::fill( lPadHisto, this->grid->nDim, this->grid->neighborhoodSize );
		//ITL_util<int>::fill( hPadHisto, this->grid->nDim, this->grid->neighborhoodSize );

		// Initialize the padded scalar field for histogram bins
		if( this->binData == NULL )
			this->binData = new ITL_field_regular<int>( this->dataField1->grid->nDim,
												this->dataField1->grid->low, this->dataField1->grid->high,
												this->dataField1->grid->lowPad, this->dataField1->grid->highPad,
												//lPadHisto, hPadHisto,
												this->dataField1->grid->neighborhoodSize );
		
		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		for( int z=0; z<this->dataField1->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField1->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField1->grid->dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					*nextVField1 = this->dataField1->datastore->array[index1d];
					*nextVField2 = this->dataField2->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location from both fields individually
					int binValue1 = histogram->get_bin_number_3D( *nextVField1, nBin );
					int binValue2 = histogram->get_bin_number_3D( *nextVField2, nBin );

					// Combine the two bin indices from the two fields
					int combinedBin = binValue2 * nBin + binValue1;

					// Store the bin ID
					this->binData->setDataAt( index1d, combinedBin );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

		delete nextVField1;
		delete nextVField2;
	   // delete lPadHisto;
	   // delete hPadHisto;

	}// end function

	void computeJointHistogramFrequencies( int nBin )
	{
		// Scan through bin Ids and keep count
		if( jointFreqList == NULL ) jointFreqList = new int[nBin*nBin];
		for( int i=0; i<nBin*nBin; i++ )
			jointFreqList[i] = 0;
		for( int i=0; i<binData->grid->nVertices; i++ )
			jointFreqList[ binData->datastore->array[i] ] ++;
	}

	/**
	 * Global joint Entropy computation function.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBins Number of bins used in histogram computation.
	 */
	void computeGlobalJointEntropyOfField( int nBin, bool toNormalize, int method = 0 )
	{
		if( method == 0 )
		{
			// Check if frequencies are already computed from bin Ids
			// If not, only then compute frequencies from bin Ids
			computeJointHistogramFrequencies( nBin*nBin );

			// Compute entropy from joint frequency list
			this->globalJointEntropy = ITL_entropycore::computeEntropy_HistogramBased( this->binData->datastore->array,
																					   this->jointFreqList,
																					   this->binData->grid->nVertices,
																					   nBin*nBin, toNormalize );

		}

	}// end function

	/**
	  * Joint historgam range set function
	  * Sets the range over which historgam computation is to be done.
	  * @param minR1 lower limit of the range along first dimension.
	  * @param maxR1 upper limit of the range along first dimension.
	  * @param minR2 lower limit of the range along second dimension.
	  * @param maxR2 upper limit of the range along second dimension.
	  */
	void setJointHistogramRange( T minR1, T maxR1, T minR2, T maxR2 )
	{
		histogramMin1 = minR1;
		histogramMax1 = maxR1;
		histogramMin2 = minR2;
		histogramMax2 = maxR2;
		histogramRangeSet = true;
	}// end function

	/**
	 * Joint histogram frequency data accessor function.
	 * Returns pointer to integer array storing the joint histogram freqencies.
	 */
	void getHistogramFrequencies( int nBin, int *jointfreqlist )
	{
		assert( this->jointFreqList != NULL );

		// Allocate memorry for the bin frequencies
		if( jointfreqlist == NULL ) jointfreqlist = new int[nBin];

		// Copy frequency
		memcpy( jointfreqlist, this->jointFreqList, sizeof( int ) * nBin );	

	}// end function

	/**
	 * Entropy accessor function.
	 * @return computed entropy.
	 */
	float getGlobalJointEntropy()
	{
		return this->globalJointEntropy;
	}// end function

};

#endif
/* ITL_ENTROPY_H_ */
