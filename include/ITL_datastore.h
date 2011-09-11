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
public:

	FILE* datafile; /**< Pointer to data file. Required if data has to be read directly from file. */
	int nel;		/**< Total number of elements in the 1D array. */
	T* array;		/**< Pointer to the 1D array containing data elements. */

public:

	/**
	 * Default Constructor.
	 */
	ITL_datastore()
	{
		this->array = NULL;
		this->datafile = NULL;

	}// end constructor 1

	/**
	 * Constructor.
	 * @param data pointer to template array conaining data.
	 */
	ITL_datastore( T* data, int nData )
	{
		//this->array = data;
		this->array = new T[nData];
		memcpy( this->array, data, sizeof(T)*nData );
		this->datafile = NULL;

	}// end constructor 1

	/**
	 * Constructor.
	 * @param nData Number of data elements.
	 * This constructor allocates memory for data.
	 */
	ITL_datastore( int nData )
	{
		this->array = new T[nData];
		this->datafile = NULL;

	}// end constructor 1

	/**
	 * Boolean function that returns TRUE if datastore contains data.
	 * @return Boolean flag.
	 */
	virtual bool hasData()
	{
		return array != NULL;
	}// end function

	/**
	 * Data accessor function.
	 * @return Pointer to data.
	 */
	virtual T* getData()
	{
		return array;
	}

	/**
	 * Destructor
	 */
	~ITL_datastore()
	{
		if( this->array != NULL ) delete [] this->array;
	}// end destructor

};

#endif
/* ITL_DATASTORE_H_ */
