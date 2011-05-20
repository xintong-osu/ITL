/**
 * Entropy computation class.
 * Created on: Nov 19, 2010.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_ENTROPY_H_
#define ITL_ENTROPY_H_

#include "ITL_header.h"
#include "ITL_field_regular.h"

template <class T>
class ITL_entropy
{
public:

	ITL_field_regular<T> *dataField;
	ITL_field_regular<int>* binData;	/**< A scalar field containing histogram bins corresponding to field points. */
	ITL_field_regular<float> *entropyField;

	float entropy;		/**< Value of computed entropy at a specified point in the field. */
	float* probarray;   /**< Probability array used in emtropy computation. */
public:

	/**
	 * Default Constructor.
	 */
	ITL_entropy( ITL_field_regular<T> *f )
	{
		this->dataField = f;
		this->binData = NULL;
		this->entropyField = NULL;
	}

	/**
	 * Histogram bin assignment function.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void computeHistogramBinField( int nBins )
	{
		assert( this->dataField->datastore->array != NULL );
	    T* nextV = new T();

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
					this->binData->setDataAt( index1d, ITL_histogram::get_bin_number_3D( (VECTOR3)*nextV, nBins ) );

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
	void computeEntropyOfField( int nBins )
	{
	    // Allocate memory for entropy field (a non-padded scalar field), if not already done
	    if( this->entropyField == NULL )
	    	this->entropyField = new ITL_field_regular<float>( this->dataField->grid->nDim,
	    													   this->dataField->grid->low, this->dataField->grid->high );

		// Compute total number of vertices in the neighborhood, including the vertext itself
		//int nNeighbors = (int)std::pow( 2.0f*this->dataField->grid->neighborhoodSize + 1.0f, this->dataField->grid->nDim );
	    int nNeighbors = (int)( ( 2.0f*this->dataField->grid->neighborhoodSizeArray[0] + 1.0f) *
	    						( 2.0f*this->dataField->grid->neighborhoodSizeArray[1] + 1.0f) *
	    						( 2.0f*this->dataField->grid->neighborhoodSizeArray[2] + 1.0f) );

	    // Compute histogram if it is not already done
		if( this->binData == NULL )
				this->computeHistogramBinField( nBins );

		int index1d = 0;
		for( int z=0; z<this->entropyField->grid->dim[2]; z++ )
		{
			for( int y=0; y<this->entropyField->grid->dim[1]; y++ )
			{
				for( int x=0; x<this->entropyField->grid->dim[0]; x++ )
				{
					// Compute and store the value of entropy at this point
					this->computeEntropySinglePoint( x, y, z, nNeighbors, nBins, index1d );

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
	void computeEntropySinglePoint( int x, int y, int z, int nNeighbors, int nBins, int entropyFieldIndex )
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
		float entropy = this->computeEntropy( binArray, nNeighbors, nBins );

		// Store entropy
		this->entropyField->setDataAt( entropyFieldIndex, entropy );

		delete binArray;

	}// end function

	/**
	 * Entropy computation function.
	 * @param binIds Pointer to array containing histogram bin assignments of the points in the neighborhood.
	 * @param nels Number of points in the neighborhood. Same as the length of bin array.
	 * @param nbins Number of bins used in histogram computation.
	 */
	float computeEntropy( int* binIds, int nels, int nbins )
	{
		// Initialize probability array
		float* probarray = new float[nbins];
		for( int i=0; i<nbins; i++ )
			probarray[i] = 0;

		// Scan through bin Ids and keep count
		for( int i = 0; i<nels; i++ )
			probarray[ binIds[i] ] ++;

		// Turn count into probabilities
		for( int i = 0; i<nbins; i++ )
			probarray[i] /= (float)nels;


		// Compute entropy
		float entropy = 0;
		for( int i = 0; i<nbins; i++ )
		{
			entropy += ( probarray[i] * ( probarray[i] == 0 ? 0 : ( log( probarray[i] ) / log(2.0) ) ) );
		}

		entropy = -entropy;

		delete [] probarray;

		return entropy;
	}
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
/* ITL_ENTROPY_H_ */
