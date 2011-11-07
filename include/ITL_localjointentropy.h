/**
 * Joint entropy computation class.
 * Created on: Jan 05, 2011.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_LOCALJOINTENTROPY_H_
#define ITL_LOCALJOINTENTROPY_H_

#include "ITL_header.h"
#include "ITL_field_regular.h"
#include "ITL_entropycore.h"

template <class T>
class ITL_localjointentropy
{
public:

	ITL_field_regular<T> *dataField1;
	ITL_field_regular<T> *dataField2;
	ITL_field_regular<int> *binData;
	ITL_field_regular<float> *jointEntropyField;

	T histogramMin1, histogramMax1, histogramMin2, histogramMax2;
	bool histogramRangeSet;

	ITL_histogram *histogram;	

public:

	/**
	 * Default Constructor.
	 */
	ITL_localjointentropy( ITL_field_regular<T> *f1, ITL_field_regular<T> *f2, ITL_histogram* hist )
	{
		this->dataField1 = f1;
		this->dataField2 = f2;
		this->binData = NULL;
		this->jointEntropyField = NULL;
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
		cout << "Scalar" << endl;
		assert( this->dataField1->datastore->array != NULL );
		assert( this->dataField2->datastore->array != NULL );
		T nextVField1, nextVField2;
		T minValue1, maxValue1, minValue2, maxValue2;

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
					nextVField2 = this->dataField1->datastore->array[index1d];

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

	/**
	 * Local joint Entropy computation function.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBins Number of bins used in histogram computation.
	 */
	void computeLocalJointEntropyOfField( int nBins )
	{
	    // Allocate memory for entropy field (a non-padded scalar field), if not already done
	    if( this->jointEntropyField == NULL )
	    	this->jointEntropyField = new ITL_field_regular<float>( this->dataField1->grid->nDim,
	    													   this->dataField1->grid->low, this->dataField1->grid->high );

		// Compute total number of vertices in the neighborhood, including the vertext itself
		int nNeighbors = (int)std::pow( 2.0f*this->dataField1->grid->neighborhoodSize + 1.0f, this->dataField1->grid->nDim );

		int index1d = 0;
		for( int z=0; z<this->jointEntropyField->grid->dim[2]; z++ )
		{
			for( int y=0; y<this->jointEntropyField->grid->dim[1]; y++ )
			{
				for( int x=0; x<this->jointEntropyField->grid->dim[0]; x++ )
				{
					// Compute and store the value of entropy at this point
					this->computeJointEntropySinglePoint( x, y, z, nNeighbors, nBins, index1d );

					// increment to the next element
					index1d ++;
				}
			}
		}

	}// end function

	/**
	 * Entropy computation function.
	 * Computes entropy at a spatial point in the field.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 * @param nNeighbors total number of neighbors in the neighborhood (including self).
	 * @param nBins Number of bins to use in histogram computation.
	 * @param entropyFieldIndex 1D index to the entopy field.
	 */
	void computeJointEntropySinglePoint( int x, int y, int z, int nNeighbors, int nBins, int entropyFieldIndex )
	{
		int* binArray = new int[nNeighbors];
		int index1d = 0;
		int binArrayIndex = 0;

		for( int k = -this->binData->grid->neighborhoodSize; k <= this->binData->grid->neighborhoodSize; k++ )
		{
			for( int j = -this->binData->grid->neighborhoodSize; j <= this->binData->grid->neighborhoodSize; j++ )
			{
				for( int i = -this->binData->grid->neighborhoodSize; i <= this->binData->grid->neighborhoodSize; i++ )
				{
					// Convert the 1D index corresponding to the point
					#ifdef DUPLICATE
					index1d = this->binData->grid->convert3DIndex( ITL_util<int>::clamp( x+i+this->binData->grid->lowPad[0], 0, this->binData->grid->dimWithPad[0]-1 ),
													ITL_util<int>::clamp( y+j+this->binData->grid->lowPad[1], 0, this->binData->grid->dimWithPad[1]-1 ),
													ITL_util<int>::clamp( z+k+this->binData->grid->lowPad[2], 0, this->binData->grid->dimWithPad[2]-1 ) );
					#endif

					#ifdef MIRROR
					index1d = this->binData->grid->convert3DIndex( ITL_util<int>::mirror( x+i+this->binData->grid->lowPad[0], 0, this->binData->grid->dimWithPad[0]-1 ),
													ITL_util<int>::mirror( y+j+this->binData->grid->lowPad[1], 0, this->binData->grid->dimWithPad[1]-1 ),
													ITL_util<int>::mirror( z+k+this->binData->grid->lowPad[2], 0, this->binData->grid->dimWithPad[2]-1 ) );
					#endif

					// Obtain the binID corresponding to the value at this location
					binArray[binArrayIndex] = this->binData->getDataAt( index1d );

					// Move to next point in the neighborhood
					binArrayIndex ++;
				}
			}
		}

		// Compute entropy
		int* localJointFreqList = NULL;
		float entropy = this->computeJointEntropy( binArray, localJointFreqList, nNeighbors, nBins*nBins );

		// Store entropy
		this->jointEntropyField->setDataAt( entropyFieldIndex, entropy );

		delete [] localJointFreqList;
		delete binArray;

	}// end function

	/**
	 * Joint entropy computation function.
	 * @param binIds Pointer to array containing histogram bin assignments of the points in the neighborhood.
	 * @param nels Number of points in the neighborhood. Same as the length of bin array.
	 * @param nbins Number of bins used in histogram computation.
	 */
	float computeJointEntropy( int* binIds, int* jointFreqArray, int nels, int nbins )
	{
		// Initialize probability array
		if( jointFreqArray == NULL ) jointFreqArray = new int[nbins];
		float* probarray = new float[nbins];
		for( int i=0; i<nbins; i++ )
		{
			jointFreqArray[i] = 0;
			probarray[i] = 0;
		}

		// Scan through bin Ids and keep count
		for( int i = 0; i<nels; i++ )
			jointFreqArray[ binIds[i] ] ++;

		// Turn count into probabilities
		for( int i = 0; i<nbins; i++ )
			probarray[i] = jointFreqArray[i] / (float)nels;


		// Compute entropy
		float entropy = 0;
		for( int i = 0; i<nbins; i++ )
		{
			entropy += ( probarray[i] * ( probarray[i] == 0 ? 0 : ( log( probarray[i] ) / log(2.0) ) ) );
		}

		entropy = -entropy;

		delete [] jointFreqArray;
		delete [] probarray;

		return entropy;
	}

	/**
	 * Entropy field accessor function.
	 * Returns computed entropy field.
	 * @return pointer to entropy field.
	 */
	ITL_field_regular<float>* getLocalJointEntropyField()
	{
		return this->jointEntropyField;
	}// end function
};

#endif
/* ITL_ENTROPY_H_ */
