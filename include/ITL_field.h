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

	ITL_grid<T> grid;			/**< Grid associated to the field. */
	ITL_datastore<T> datastore;		/**< Data store associated to the field. */

	//ITL_grid<T> *grid;			/**< Grid associated to the field. */
	//ITL_datastore<T> *datastore;		/**< Data store associated to the field. */
	//ITL_field<T>* subfieldArray;		/**< An array of sub fields (also referred to as blocks). This can be potentially useful for recursive subdivision of field. */
	//int* nPartition;			/**< Number of subfields along each dimension. */

public:

	/**
	 * Default constructor.
	 */
	ITL_field (){}

	/**
	 * Destructor.
	 */
	virtual ~ITL_field()
	{
	}

	/**
	 * Pure virtual function for regular partitioning of the field.
	 */
	virtual void
	partitionField( int* nblocks )
	{
	}
	/**
	 * Pure virtual function for creating a partition or subfield
	 * within the specified bounds.
	 */
	virtual ITL_field*
	createSubField( float* low, float* high )
	{
		return NULL;
	}


	/**
	 * Pure virtual function for setting the bounds of a field with no ghost layers.
	 * Calls corresponding function of the grid.
	 */
	virtual void
	setBounds( float*, float* )
	{
	}

	/**
	 * Pure virtual function for setting the bounds of a grid with ghost layers.
	 * Calls corresponding function of the grid.
	 */
	virtual void
	setBounds( float*, float*, int*, int*, int )
	{
	}

	/**
	 * Data mutator function type 1.
	 */
	virtual void
	setDataAt( int id, T data )
	{
	}

	/**
	 * Data mutator function type 2.
	 */
	virtual void
	setDataAt( int x, int y, int z, T data )
	{
	}

	/**
	 * Data mutator function type 3.
	 */
	virtual void
	setDataBetween( T* data )
	{
	}

	/**
	 * Data mutator function type 4.
	 */
	virtual void setDataFull( T* data ){}

	/**
	 * Field bounds accessor function type 1.
	 */
	virtual void
	getBounds( int *limits )
	{
	}

	/**
	 * Field bounds accessor function type 2.
	 */
	virtual void
	getBounds( int* l, int* h )
	{
	}

	/**
	 * Field bounds accessor function type 3.
	 */
	virtual void
	getBounds( float* l, float* h )
	{
	}

	/**
	 * Field size (in terms of total number of grid vertices) accessor function.
	 */
	virtual int
	getSize()
	{
		return 0;
	}

	/**
	 * Field size (in terms of number of points along each dimension) accessor function.
	 */
	virtual void
	getSize( int* d )
	{
	}

	/**
	 * Field size (in terms of number of points with ghost layers along each dimension) accessor function.
	 */
	virtual void
	getSizeWithPad( int* d )
	{
	}

	/**
	 * Pad size accessor function.
	 */
	virtual void
	getPadSize( int* lowpad, int* highpad )
	{
	}

	/**
	 * Neighborhood size accessor function.
	 */
	virtual int
	getNeighborhoodSize()
	{
		return 0;
	}

	/**
	 * Neighborhood size array accessor function.
	 */
	virtual void
	getNeighborhoodSize( int* neighborsize )
	{
	}

	/**
	 * Number of dimensions (max 4) accessor function.
	 */
	virtual int
	getNumDim()
	{
		return 0;
	}

	/**
	 * Data accessor function type 1.
	 */
	virtual T
	getDataAt( int id )
	{
		return (T)0.0f;
	}

	/**
	 * Data accessor function type 2.
	 */
	virtual T
	getDataAt( int x, int y, int z )
	{
		return (T)0.0f;
	}

	/**
	 * Data accessor function type 3.
	 */
	virtual void
	getDataBetween( int* lowBoundary, int* highBoundary, T* retData )
	{
	}

	/**
	 * Data accessor function type 4.
	 */
	virtual void getDataBetween( float* lowBoundary, float* highBoundary, T* retData )
	{
	}

	/**
	 * Data accessor function type 5.
	 */
	virtual T*
	getDataFull()
	{
		return NULL;
	}

	/**
	 * Grid accessor function.
	 */
	virtual ITL_grid<T>*
	getGrid()
	{
		return &grid;
	}

};

#endif
/* ITL_FIELD_H_ */
