/**
 * The container class for data.
 * This template class is the container for the raw data for a field. This class does not have any
 * information regarding indexing and spatial arrangement of the data.
 * Created on: Nov 18, 2010.
 * @authors Abon
 * @author Teng-Yok
 * @see ITL_grid
 */

#ifndef ITL_DATASTORE_H_
#define ITL_DATASTORE_H_

#include "ITL_header.h"

template <class T>
class ITL_datastore
{
	//FILE* datafile; /**< Pointer to data file. Required if data has to be read directly from file. */
	int nel;		/**< Total number of elements in the 1D array. */
	T* array;		/**< Pointer to the 1D array containing data elements. */

public:

	/**
	 * Default Constructor.
	 */
	ITL_datastore()
	{
		this->array = NULL;

	}// end constructor 1

	/**
	 * Constructor.
	 * @param data pointer to template array conaining data.
	 */
	ITL_datastore( T* data, int ndata )
	{
		//this->array = data;
		nel = ndata;
		this->array = NULL;
		this->array = new T[nel];
		memcpy( this->array, data, sizeof(T)*nel );

	}// end constructor 1

	ITL_datastore( const ITL_datastore<T>& that )
	{
		this->nel = that.nel;
		this->array = NULL;
		this->array = new T[nel];
		memcpy( this->array, that.array, sizeof(T)*nel );
	}

	ITL_datastore( ITL_datastore<T>& that )
	{
		this->nel = that.nel;
		this->array = NULL;
		this->array = new T[nel];
		memcpy( this->array, that.getData(), sizeof(T)*nel );
	}

	ITL_datastore& operator= ( const ITL_datastore<T>& that )
	{
		if ( this != &that ) // protect against invalid self-assignment
		{
			this->nel = that.nel;
			this->array = NULL;
			this->array = new T[nel];
			memcpy( this->array, that.array, sizeof(T)*nel );
		}
		// by convention, always return *this
		return *this;
	}

	ITL_datastore& operator= ( ITL_datastore<T>& that )
	{
		if ( this != &that ) // protect against invalid self-assignment
		{
			this->nel = that.nel;
			this->array = NULL;
			this->array = new T[nel];
			memcpy( this->array, that.array, sizeof(T)*nel );
		}
		// by convention, always return *this
		return *this;
	}


	/**
	 * Destructor
	 */
	virtual
	~ITL_datastore()
	{
		if( this->array != NULL ) delete [] this->array;
	}// end destructor

	void
	init( int ndata )
	{
		nel = ndata;
		this->array = new T[nel];

	}// end constructor 1

	void
	init( T* data, int ndata )
	{
		nel = ndata;
		this->array = new T[nel];
		memcpy( this->array, data, sizeof(T)*nel );

	}// end constructor 1

	/**
	 * Boolean function that returns TRUE if datastore contains data.
	 * @return Boolean flag.
	 */
	virtual bool hasData()
	{
		return array != NULL;
	}// end function

	virtual void
	setDataAt( T v, int i )
	{
		array[i] = v;
	}

	virtual void
	setDataFull( T* ptr )
	{
		array = ptr;
	}

	virtual void
	setDataFull( T* ptr, int ndata )
	{
		nel = ndata;
		if( this->array != NULL )
			delete [] this->array;
		this->array = new T[nel];
		memcpy( this->array, ptr, sizeof(T)*nel );
	}


	/**
	 * Data accessor function.
	 * @return Pointer to data.
	 */
	virtual T*
	getData()
	{
		return array;
	}

	virtual int
	getDataSize()
	{
		return nel;
	}

	virtual T
	getDataAt( int i )
	{
		assert( array != NULL );
		assert( i < nel );
		return array[i];
	}

};

#endif
/* ITL_DATASTORE_H_ */
