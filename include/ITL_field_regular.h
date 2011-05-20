/**
 * Regular field class inherited from ITL_field.
 * Container for regularly arranged static/time-varying scalar/vector data.
 * Created on: Dec 09, 2010
 * @author Abon
 * @author Teng-Yok
 * @see ITL_field
 */

#ifndef ITL_FIELD_REGULAR_H_
#define ITL_FIELD_REGULAR_H_

#include "ITL_field.h"
#include "ITL_grid_regular.h"

template <class T>
class ITL_field_regular: public ITL_field<T>
{
public:

	/**
	 * Default constructor.
	 */
	ITL_field_regular()
	{
	}

	/**
	 * Constructor type 2.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 */
	ITL_field_regular( T* data, int ndim, int* dim )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds: 0 to N-1
		float* lowEnd = new float[ndim];
		float* highEnd = new float[ndim];
		for( int i=0; i<ndim; i++ )
		{
			lowEnd[i] = 0.0f;
			highEnd[i] = (float)(dim[i]-1);
		}
		this->setBounds( lowEnd, highEnd );

		// Initialize data store
		this->datastore = new ITL_datastore<T>( data );

	}// constructor 1

	/**
	 * Constructor type 3.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	ITL_field_regular( T* data, int ndim, float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsize );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( data );

	}// Constructor

	ITL_field_regular( T* data, int ndim, float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( data );

	}// Constructor

	/**
	 * Constructor type 3.
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	ITL_field_regular( int ndim, float* l, float* h )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( this->grid->dimWithPad, this->grid->nDim ) );

	}// constructor

	/**
	 * Constructor type 4.
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	ITL_field_regular( int ndim, float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsize );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( this->grid->dimWithPad, this->grid->nDim ) );

	}// Constructor

	ITL_field_regular( int ndim, float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( this->grid->dimWithPad, this->grid->nDim ) );

	}// Constructor


	/**
	 * Pure virtual function for setting the bounds of a field with no ghost layers.
	 * Calls corresponding function of the grid.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	void setBounds( float* l, float* h )
	{
		// call same method of the grid
		this->grid->setBounds( l, h );

	}// end function

	/**
	 * Pure virtual function for setting the bounds of a field with ghost layers.
	 * Calls corresponding function of the grid.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	void setBounds( float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// call same method of the grid
		this->grid->setBounds( l, h , lPad, hPad, neighborhoodsize );

	}// end function

	void setBounds( float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// call same method of the grid
		this->grid->setBounds( l, h , lPad, hPad, neighborhoodsizearray );

	}// end function

	/**
	 * Function for creating a partition or subfield
	 * within the specified bounds.
	 * @param nblocks Number of partitions along each dimension.
	 */
	virtual void partitionField ( int* nblocks )
	{
		this->nPartition = nblocks;
		this->blockSize = new float[this->grid->nDim];
		float* lowSub = new float[this->grid->nDim];
		float* highSub = new float[this->grid->nDim];

		// Compute the dimension of each block
		this->blockSize[0] = (this->grid->dim[0]-1) / (float)this->nPartition[0];
		this->blockSize[1] = (this->grid->dim[1]-1) / (float)this->nPartition[1];
		this->blockSize[2] = (this->grid->dim[2]-1) / (float)this->nPartition[2];

		#ifdef DEBUG_MODE
		printf( "Number of partitions (ie, subfields) to be created: %d, %d, %d\n", this->nPartition[0], this->nPartition[1], this->nPartition[2] );
		printf( "Dimension of a subfield: %f, %f, %f\n", this->blockSize[0], this->blockSize[1], this->blockSize[2] );
		#endif

		// Compute total number of blocks
		int nTotalBlocks = ITL_util<int>::prod( this->nPartition, this->grid->nDim );

		// Allocate memory for array of subfields
		this->subfieldArray = new ITL_field_regular<T>[nTotalBlocks];

		int blockIndex = 0;

		for( int z=0; z<this->nPartition[2]; z++ )
		{
			lowSub[2] = z * this->blockSize[2];
			highSub[2] = lowSub[2] + this->blockSize[2];
			for(int y=0; y<this->nPartition[1]; y++ )
			{
				lowSub[1] = y * this->blockSize[1];
				highSub[1] = lowSub[1] + this->blockSize[1];
				for( int x=0; x<this->nPartition[0]; x++)
				{
					lowSub[0] = x * this->blockSize[0];
					highSub[0] = lowSub[0] + this->blockSize[0];

					// Initialize sub field
					this->subfieldArray[blockIndex].setBounds( lowSub, highSub,
															   this->grid->lowPad, this->grid->highPad,
															   this->grid->neighborhoodSize );

					// load Data to the subfield
					this->subfieldArray[blockIndex].setDataFull( this->getDataBetween( lowSub, highSub ) );

					// Increment to the next block
					blockIndex ++;
				}
			}
		}

		delete [] lowSub;
		delete [] highSub;
	}// end function

	/**
	 * Function for creating a partition or subfield
	 * within the specified bounds.
	 * @param low Lower bound along each dimension.
	 * @param high Higher bound along each dimension.
	 */
	ITL_field<T>* createSubField( float* low, float* high )
	{
		// Compute the dimension of each block
		this->blockSize[0] = (high[0]-low[0]);
		this->blockSize[1] = (high[1]-low[1]);
		this->blockSize[2] = (high[2]-low[2]);

		#ifdef DEBUG_MODE
		printf( "Dimension of newly created subfield: %f, %f, %f\n", this->blockSize[0], this->blockSize[1], this->blockSize[2] );
		#endif

		// Allocate memory for 3D block
		ITL_field_regular<T>* newField = new ITL_field_regular<T>( this->getDataBetween( low, high ), this->grid->nDim, this->grid->dim );

		return newField;

	}// end function

	/**
	 * Data accessor function type 1.
	 * @param id 1D index of data array
	 * @return field value.
	 */
	virtual T getDataAt( int id )
	{
		// Get data
		return this->datastore->array[id];
	}// end function

	/**
	 * Data accessor function type 2.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 * @return field value.
	 */
	virtual T getDataAt( int x, int y, int z )
	{
		// Convert index to 1D
		int index1d = this->grid->convert3DIndex( x, y, z );

		// Get data
		return this->datastore->array[index1d];
	}// end function

	/**
	 * Data accessor function type 3.
	 * Returns chunk of data within specified bound.
	 * @param lowBoundary Lower bound along each dimension.
	 * @param highBoundary Higher bound along each dimension.
	 * @return pointer to field value.
	 */
	virtual T* getDataBetween( float* lowBoundary, float* highBoundary  )
	{
		int* lowBoundaryInt = new int[this->grid->nDim];
		int* highBoundaryInt = new int[this->grid->nDim];

		for( int i=0; i<this->grid->nDim; i++ )
		{
			lowBoundaryInt[i] = (int)floor( lowBoundary[i] );
			highBoundaryInt[i] = (int)ceil( highBoundary[i] );
		}

		return this->getDataBetween( lowBoundaryInt, highBoundaryInt );

	}// end function

	/**
	 * Data accessor function type 4.
	 * Returns chunk of data within specified bound.
	 * @param lowBoundary Lower bound along each dimension.
	 * @param highBoundary Higher bound along each dimension.
	 * @return pointer to field value.
	 */
	virtual T* getDataBetween( int* lowBoundary, int* highBoundary  )
	{
		// Count number of vertices requested
		int* dimLength = ITL_util<int>::subtractArrays( highBoundary, lowBoundary, this->grid->nDim );
		ITL_util<int>::addArrayScalar( dimLength, 1, this->grid->nDim );
		int nV = ITL_util<int>::prod( dimLength, this->grid->nDim );

		// Allocate memory to return data
		T* retData = new T[nV];

		int globalOffset = 0;
		int localOffset = 0;
		for( int z=0; z<dimLength[2]; z++ )
		{
			for( int y=0; y<dimLength[1]; y++ )
			{
				globalOffset = this->grid->convert3DIndex( z+lowBoundary[2], y+lowBoundary[1], lowBoundary[0] );

				memcpy( (T*)(retData + localOffset) , (T*)(this->datastore->array + globalOffset) , dimLength[0] * sizeof( T ) );
				localOffset += dimLength[0];
			}
		}

		return retData;

	}// end function

	/**
	 * Data accessor function type 5.
	 * Returns entire field data.
	 * @return pointer to field value.
	 */
	virtual T* getDataFull()
	{
		return this->datastore->array;
	}// end function

	/**
	 * Data mutator function type 1.
	 * Sets data at a particular array location.
	 * @param id 1D index of data array
	 * @param data field value.
	 */
	virtual void setDataAt( int id, T data )
	{
		// Set data
		this->datastore->array[id] = data;
	}// end function

	/**
	 * Data mutator function type 2.
	 */
	virtual void setDataAt( int x, int y, int z, T data )
	{
		// Convert index to 1D
		int index1d = this->grid->convert3DIndex( x, y, z );

		// Set data
		this->datastore->array[index1d] = data;
	}// end function

	/**
	 * Data mutator function type 3.
	 */
	virtual void setDataBetween( T* data )
	{
		// Not implemented yet
	}// end function

	/**
	 * Data mutator function type 4.
	 * Sets data for entire field.
	 * @param data pointer to data
	 */
	virtual void setDataFull( T* data )
	{
		this->datastore->array = data;
	}// end function

	/**
	 * Destructor.
	 */
	~ITL_field_regular()
	{
	}


};

#endif
/* ITL_FIELD_REGULAR_H_ */
