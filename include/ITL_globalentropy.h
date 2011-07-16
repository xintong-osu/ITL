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

	ITL_field_regular<T> *dataField;
	ITL_field_regular<int>* binData;		/**< A scalar field containing histogram bins corresponding to field points. */
	int *freqList;

	T histogramMin;
	T histogramMax;
	bool histogramRangeSet;

	float globalEntropy;				/**< Value of computed global entropy of the field. */
	float* probarray;   				/**< Probability array used in entropy computation. */

public:

	/**
	 * Default Constructor.
	 */
	ITL_globalentropy( ITL_field_regular<T> *f )
	{
		dataField = f;
		binData = NULL;
		freqList = NULL;
		histogramRangeSet = false;

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
			
			computeHistogramBinField_Vector( nBin, binMapFile );
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
			//this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim,
			//											this->dataField->grid->low, this->dataField->grid->high,
			//											this->dataField->grid->lowPad, this->dataField->grid->highPad,
			//											//lPadHisto, hPadHisto,
			//											this->dataField->grid->neighborhoodSize );
			this->binData = new ITL_field_regular<int>( this->dataField->grid->nDim, this->dataField->grid->low, this->dataField->grid->high, this->dataField->grid->lowPad, this->dataField->grid->highPad, this->dataField->grid->neighborhoodSizeArray );
			//lPadHisto, hPadHisto,
											
		
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
		printf( "Min: %f Max: %f Range: %f of the scalar values\n", minValue, maxValue, rangeValue );
		printf( "Binwidth: %f\n", binWidth );		
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
	void computeHistogramBinField_Vector( int nBin, char* binMapFile = NULL )
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
	 * @param nBins Number of bins used in histogram computation/Number of sample points in KDE estimation.
	 * @paran toNormalize boolean flag, TRUE indicates that the computed entropy value needs to be normalized
	 */
	void computeGlobalEntropyOfField( int nBin, bool toNormalize, int method = 0 )
	{
		if( method == 0 )
			globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( this->binData->datastore->array, this->freqList, this->binData->grid->nVertices, nBin, toNormalize );		
		else
			//globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( this->binData->datastore->array, this->binData->grid->nVertices, nBin, toNormalize );		
			globalEntropy = ITL_entropycore::computeEntropy_KDEBased( this->dataField->datastore->array, this->dataField->grid->nVertices, 0.0, toNormalize );
	
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

		// Allocate memorry for the bin frequencies
		if( freqlist == NULL ) freqlist = new int[nBin];

		// Copy frequency
		memcpy( freqlist, this->freqList, sizeof( int ) * nBin );
		
	

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
