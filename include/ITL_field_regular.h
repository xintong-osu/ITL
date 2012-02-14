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

#include "ITL_util.h"
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
		this->grid = NULL;
		this->datastore = NULL;
	}

	/**
	 * Constructor. 
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
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// constructor

	/**
	 * Constructor.
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
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// Constructor

	/**
	 * Constructor.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsizearray Neighborhood length for each dimension.
	 */
	ITL_field_regular( int ndim, float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// Constructor

	/**
	 * Constructor.
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
		this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// constructor 1

	/**
	 * Constructor. 
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	ITL_field_regular( T* data, int ndim, float* l, float* h )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// constructor

	/**
	 * Constructor.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	ITL_field_regular( T* data, int ndim,
					   float* l, float* h,
					   int* lPad, int* hPad,
					   int neighborhoodsize )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsize );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// Constructor

	/**
	 * Constructor.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsizearray Neighborhood length for each dimension.
	 */
	ITL_field_regular( T* data, int ndim,
					   float* l, float* h,
					   int* lPad, int* hPad,
					   int* neighborhoodsizearray )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );

	}// Constructor

	void initialize( int ndim, float *l, float *h )
	{
		// Initialize grid
		this->grid = new ITL_grid_regular<T>( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dim[0], ndim ) );	
	}


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
	 * @param nblocks Number of blocks along each dimension.
	 * @paam subfieldArray A null pointer or a pointer to sub-blocks to be created.  
	 * @param isSharingNeighbor Flag that denotes if two adjacent blocks are sharing a common neighbor. Default value is false.
	 */
	virtual void partitionField ( int* nBlock, ITL_field_regular<T> *subfieldArray, bool isSharingNeighbor = true )
	{		
		float *blockSize = new float[this->grid->nDim];
		float* lowSub = new float[this->grid->nDim];
		float* highSub = new float[this->grid->nDim];

		// Compute the dimension of each block
		blockSize[0] = ceil( (this->grid->dim[0]-1) / (float)nBlock[0] );
		blockSize[1] = ceil( (this->grid->dim[1]-1) / (float)nBlock[1] );
		blockSize[2] = ceil( (this->grid->dim[2]-1) / (float)nBlock[2] );
		//blockSize[0] = floor( (this->grid->dim[0]-1) / (float)nBlock[0] );
		//blockSize[1] = floor( (this->grid->dim[1]-1) / (float)nBlock[1] );
		//blockSize[2] = floor( (this->grid->dim[2]-1) / (float)nBlock[2] );


		#ifdef DEBUG_MODE
		printf( "Dimension of field: %d, %d, %d\n", this->grid->dim[0], this->grid->dim[1], this->grid->dim[2] );
		printf( "Number of blocks (ie, subfields) to be created: %d, %d, %d\n", nBlock[0], nBlock[1], nBlock[2] );
		printf( "Dimension of a subfield: %f, %f, %f\n", blockSize[0], blockSize[1], blockSize[2] );
		#endif

		// Compute total number of blocks
		int nTotalBlocks = ITL_util<int>::prod( nBlock, this->grid->nDim );

		// Allocate memory for array of subfields
		if( subfieldArray == NULL )	subfieldArray = new ITL_field_regular<T>[nTotalBlocks];

		// Partitioning loop
		int blockIndex = 0;
		for( int z=0; z<nBlock[2]; z++ )
		{
			// Determine [zmin, zmax] for the next block
			lowSub[2] = z * blockSize[2];			
			if( isSharingNeighbor == true )		highSub[2] = std::min( (float)this->grid->highInt[2], lowSub[2] + blockSize[2] );
						
			for(int y=0; y<nBlock[1]; y++ )
			{
				// Determine [ymin, ymax] for the next block
				lowSub[1] = y * blockSize[1];
				if( isSharingNeighbor == true ) highSub[1] = std::min( (float)this->grid->highInt[1], lowSub[1] + blockSize[1] );
				
				for( int x=0; x<nBlock[0]; x++)
				{
					// Determine [xmin, xmax] for the next block
					lowSub[0] = x * blockSize[0];
					if( isSharingNeighbor == true )	highSub[0] = std::min( (float)this->grid->highInt[0], lowSub[0] + blockSize[0] );
					
					// Initialize subfield 
					if( subfieldArray[blockIndex].grid == NULL )  		subfieldArray[blockIndex].grid = new ITL_grid_regular<T>( this->grid->nDim );
					if( subfieldArray[blockIndex].datastore == NULL )  	subfieldArray[blockIndex].datastore = new ITL_datastore<T>();
					subfieldArray[blockIndex].setBounds( lowSub, highSub, this->grid->lowPad, this->grid->highPad, this->grid->neighborhoodSize );
					#ifdef DEBUG_MODE
					printf( "%d: %f %f %f %f %f %f\n", blockIndex, lowSub[0], lowSub[1], lowSub[2], highSub[0], highSub[1], highSub[2] );
					#endif

					// load (copy from datastore of full field) Data to the subfield
					subfieldArray[blockIndex].setDataFull( this->getDataBetween( lowSub, highSub ) );
					//float m = ITL_util<SCALAR>::Min( (SCALAR*)subfieldArray[blockIndex].datastore->array, subfieldArray[blockIndex].grid->nVertices );
					//float M = ITL_util<SCALAR>::Max( (SCALAR*)subfieldArray[blockIndex].datastore->array, subfieldArray[blockIndex].grid->nVertices );
					//cout << m << " " << M << endl;

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
				globalOffset = this->grid->convert3DIndex( lowBoundary[0], lowBoundary[1]+y, lowBoundary[2]+z );

				memcpy( (T*)(retData + localOffset ) , (T*)(this->datastore->array + globalOffset ) , dimLength[0] * sizeof( T ) );
				localOffset += dimLength[0];
			}
		}
		
		return retData;

	}// end function
	
	virtual void getDataBetween2( float* lowBoundary, float* highBoundary, T* retData )
	{
		// ADD-BY-LEETEN 02/13/2012-BEGIN
		#ifdef	WIN32
		int* lowBoundaryInt = new int[this->grid->nDim];
		int* highBoundaryInt = new int[this->grid->nDim];
		#else	// #ifdef	WIN32
		// ADD-BY-LEETEN 02/13/2012-END

		int lowBoundaryInt[this->grid->nDim];
		int highBoundaryInt[this->grid->nDim];

		#endif	// #ifdef	WIN32		// ADD-BY-LEETEN 02/13/2012

		for( int i=0; i<this->grid->nDim; i++ )
		{
			lowBoundaryInt[i] = (int)floor( lowBoundary[i] );
			highBoundaryInt[i] = (int)ceil( highBoundary[i] );
		}

		this->getDataBetween2( lowBoundaryInt, highBoundaryInt, retData );
		
		// ADD-BY-LEETEN 02/13/2012-BEGIN
		#ifdef	WIN32
		delete [] lowBoundaryInt;
		delete [] highBoundaryInt;
		#endif	// #ifdef	WIN32
		// ADD-BY-LEETEN 02/13/2012-END
	}// end function

	/**
	 * Data accessor function type 4.
	 * Returns chunk of data within specified bound.
	 * @param lowBoundary Lower bound along each dimension.
	 * @param highBoundary Higher bound along each dimension.
	 * @return pointer to field value.
	 */
	virtual void getDataBetween2( int* lowBoundary, int* highBoundary, T* retData )
	{
		// Count number of vertices requested
		int* dimLength = ITL_util<int>::subtractArrays( highBoundary, lowBoundary, this->grid->nDim );
		ITL_util<int>::addArrayScalar( dimLength, 1, this->grid->nDim );
		int nV = ITL_util<int>::prod( dimLength, this->grid->nDim );

		int globalOffset = 0;
		int localOffset = 0;
		for( int z=0; z<dimLength[2]; z++ )
		{
			for( int y=0; y<dimLength[1]; y++ )
			{
				globalOffset = this->grid->convert3DIndex( lowBoundary[0], lowBoundary[1]+y, lowBoundary[2]+z );

				memcpy( (T*)(retData + localOffset ) , (T*)(this->datastore->array + globalOffset ) , dimLength[0] * sizeof( T ) );
				localOffset += dimLength[0];
			}
		}

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
		//if( this->grid != NULL )	delete this->grid;
		//if( this->datastore != NULL )	delete this->datastore;
		//delete grid;
		//delete datastore;
	}
};

#endif
/* ITL_FIELD_REGULAR_H_ */
