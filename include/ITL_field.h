/**
 * Field base class.
 * A generic class for a field which is a generic container for static/time-varying
 * scalar/vector data.
 * Created on: Nov 17, 2010.
 * @author Abon
 * @author Teng-Yok
 * @see ITL_field regular
 */

#ifndef ITL_FIELD_H_
#define ITL_FIELD_H_

#include "ITL_datastore.h"
#include "ITL_grid.h"

template <class T>
class ITL_field
{
public:

	ITL_grid<T> *grid;			/**< Grid associated to the field. */
	ITL_datastore<T> *datastore;		/**< Data store associated to the field. */
	//ITL_field<T>* subfieldArray;		/**< An array of sub fields (also referred to as blocks). This can be potentially useful for recursive subdivision of field. */
	//int* nPartition;			/**< Number of subfields along each dimension. */

public:

	/**
	 * Default constructor.
	 */
	ITL_field (){}
	/**
	 * Pure virtual function.
	 */
	virtual void setBounds(){}
	/**
	 * Pure virtual function for setting the bounds of a field with no ghost layers.
	 * Calls corresponding function of the grid.
	 */
	virtual void setBounds( float*, float* ){}
	/**
	 * Pure virtual function for setting the bounds of a grid with ghost layers.
	 * Calls corresponding function of the grid.
	 */
	virtual void setBounds( float*, float*, int*, int*, int ){}

	/**
	 * Pure virtual function for regular partitioning of the field.
	 */
	virtual void partitionField( int* nblocks ){}
	/**
	 * Pure virtual function for creating a partition or subfield
	 * within the specified bounds.
	 */
	virtual ITL_field* createSubField( float* low, float* high ){ return NULL; }
	/**
	 * Data accessor function type 1.
	 */
	virtual T getDataAt( int id ){ return (T)0.0f; }
	/**
	 * Data accessor function type 2.
	 */
	virtual T getDataAt( int x, int y, int z ){ return (T)0.0f; }
	/**
	 * Data accessor function type 3.
	 */
	virtual T* getDataBetween( int* lowBoundary, int* highBoundary  ){ return NULL; }
	/**
	 * Data accessor function type 4.
	 */
	virtual T* getDataBetween( float* lowBoundary, float* highBoundary  ){ return NULL; }
	/**
	 * Data accessor function type 5.
	 */
	virtual T* getDataFull(){ return NULL; }

	/**
	 * Data mutator function type 1.
	 */
	virtual void setDataAt( int id, T data ){}
	/**
	 * Data mutator function type 2.
	 */
	virtual void setDataAt( int x, int y, int z, T data ){}
	/**
	 * Data mutator function type 3.
	 */
	virtual void setDataBetween( T* data ){}
	/**
	 * Data mutator function type 4.
	 */
	virtual void setDataFull( T* data ){}

	/**
	 * Destructor.
	 */
	virtual ~ITL_field() {}
};

#endif
/* ITL_FIELD_H_ */
