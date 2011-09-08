/**
 * Local entropy field computation class.
 * Created on: Nov 19, 2010.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_LOCALENTROPY_H_
#define ITL_LOCALENTROPY_H_

#include "ITL_header.h"
// ADD-BY-LEETEN 07/18/2011-BEGIN
#include "ITL_histogram.h"
// ADD-BY-LEETEN 07/18/2011-END
#include "ITL_field_regular.h"
#include "ITL_entropycore.h"

template <class T>
class ITL_localentropy
{
public:

	ITL_field_regular<T> *dataField;
	ITL_field_regular<int>* binData;	/**< A scalar field containing histogram bins corresponding to field points. */
	ITL_field_regular<float> *entropyField;

	T histogramMin;
	T histogramMax;
	bool histogramRangeSet;

	float globalEntropy;				/**< Value of computed entropy at a specified point in the field. */
	float* probarray;   				/**< Probability array used in emtropy computation. */
public:

	/**
	 * Default Constructor.
	 */
	ITL_localentropy( ITL_field_regular<T> *f )
	{
		this->dataField = f;
		this->binData = NULL;
		this->entropyField = NULL;
		histogramRangeSet = false;		// ADD-BY-ABON 07/19/2011
	}// End constructor
	
	/**
	 * Histogram bin assignment function.
	 * Calls the appropriate function based on field type
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField( char *fieldType, int nBin = 0 )
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
			
			computeHistogramBinField_Vector( nBin );
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
	        T nextV, minValue, maxValue, rangeValue;
		
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

		// MOD-BY-ABON 07/19/2011-BEGIN			
		if( histogramRangeSet == false )
		{
			// Get min-max values of the scalar field
			minValue = ITL_util<T>::Min( this->dataField->datastore->array, this->dataField->grid->nVertices );
			maxValue = ITL_util<T>::Max( this->dataField->datastore->array, this->dataField->grid->nVertices );
		}
		else
		{	
			minValue = histogramMin;
			maxValue = histogramMax;
		}
		rangeValue = maxValue - minValue;

		// Compute bin width			
		float binWidth = rangeValue / (float)nBin;
		// MOD-BY-ABON 07/19/2011-END					

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
					nextV = this->dataField->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV - minValue ) / binWidth  );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					//if( binId != 0 ) cout << binId << endl;
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
	void computeHistogramBinField_Vector( int nBin )
	{
		printf( "In vector\n");
		cout << this->dataField->grid->neighborhoodSizeArray[0] << " " << this->dataField->grid->neighborhoodSizeArray[1] << " " << this->dataField->grid->neighborhoodSizeArray[2] << endl;
		
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
		for( int z=0; z<this->dataField->grid->dimWithPad[2]; z++ )
		{
			for( int y=0; y<this->dataField->grid->dimWithPad[1]; y++ )
			{
				for( int x=0; x<this->dataField->grid->dimWithPad[0]; x++ )
				{
					// Get vector at location
					*nextV = this->dataField->datastore->array[index1d];

					// Obtain the binID corresponding to the value at this location
					this->binData->setDataAt( index1d, ITL_histogram::get_bin_number_3D( *nextV, nBin ) );

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
	 * @param nBins Number of bins used in histogram computation.
	 */
	void computeEntropyOfField( int nBins, bool toNormalize )
	{
		// Allocate memory for entropy field (a non-padded scalar field), if not already done
		if( this->entropyField == NULL )
			this->entropyField = new ITL_field_regular<SCALAR>( this->dataField->grid->nDim, this->dataField->grid->low, this->dataField->grid->high );

		// Compute total number of vertices in the neighborhood, including the vertext itself
		//int nNeighbors = (int)std::pow( 2.0f*this->dataField->grid->neighborhoodSize + 1.0f, this->dataField->grid->nDim );
		int nNeighbors = (int)( ( 2.0f*this->dataField->grid->neighborhoodSizeArray[0] + 1.0f ) *
						     ( 2.0f*this->dataField->grid->neighborhoodSizeArray[1] + 1.0f ) *
						     ( 2.0f*this->dataField->grid->neighborhoodSizeArray[2] + 1.0f ) );

		// Compute histogram if it is not already done
		//if( this->binData == NULL )
		//	this->computeHistogramBinField( nBins );

		int index1d = 0;
		for( int z=0; z<this->entropyField->grid->dim[2]; z++ )
		{
			for( int y=0; y<this->entropyField->grid->dim[1]; y++ )
			{
				for( int x=0; x<this->entropyField->grid->dim[0]; x++ )
				{
					// Compute and store the value of entropy at this point
					this->computeEntropySinglePoint( x, y, z, nNeighbors, nBins, index1d, toNormalize );

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
	void computeEntropySinglePoint( int x, int y, int z, int nNeighbors, int nBins, int entropyFieldIndex, bool toNormalize )
	{
		int* binArray = new int[nNeighbors];
		int index1d = 0;
		int binArrayIndex = 0;

//		for( int k = -this->binData->grid->neighborhoodSize; k <= this->binData->grid->neighborhoodSize; k++ )
//		{
//			for( int j = -this->binData->grid->neighborhoodSize; j <= this->binData->grid->neighborhoodSize; j++ )
//			{
//				for( int i = -this->binData->grid->neighborhoodSize; i <= this->binData->grid->neighborhoodSize; i++ )
//				{
		for( int k = -this->binData->grid->neighborhoodSizeArray[2]; k <= this->binData->grid->neighborhoodSizeArray[2]; k++ )
		{
			for( int j = -this->binData->grid->neighborhoodSizeArray[1]; j <= this->binData->grid->neighborhoodSizeArray[1]; j++ )
			{
				for( int i = -this->binData->grid->neighborhoodSizeArray[0]; i <= this->binData->grid->neighborhoodSizeArray[0]; i++ )
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
		int* localFreqList = new int[nBins];
		float entropy = ITL_entropycore::computeEntropy_HistogramBased( binArray, localFreqList, nNeighbors, nBins, toNormalize );
		if( entropy > 0 ) cout << entropy << endl;

		// Store entropy
		this->entropyField->setDataAt( entropyFieldIndex, entropy );

		if( localFreqList != NULL ) delete [] localFreqList;
		delete [] binArray;

	}// end function

	// ADD-BY-ABON 07/19/2011-BEGIN	
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
	// ADD-BY-ABON 07/19/2011-END	


	/**
	 * Entropy accessor function.
	 * @return computed entropy.
	 */
	float getEntropy();

	/**
	 * Entropy field accessor function.
	 * Returns computed entropy field.
	 * @return pointer to entropy field.
	 */
	ITL_field_regular<float>* getEntropyField()
	{
		return this->entropyField;
	}// end function
};

#endif
/* ITL_LOCALENTROPY_H_ */
