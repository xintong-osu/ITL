/**
 * Grid base class.
 * A generic class for grid which contains information about spatial arrangement
 * of field data points in a Cartesian space.
 * Created on: Nov 17, 2010.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GRID_H_
#define ITL_GRID_H_

#include "ITL_util.h"

template <class T>
class ITL_grid
{
public:

	int nDim;			/**< Number of dimensions */
	int* dim;			/**< Array for length of each dimension */
	int* dimWithPad;		/**< Array for length of each dimension along with ghost layers (if any) */
	int nVertices;			/**< Number of vertices in Cartesian space. */
	int nVerticesWithPad;		/**< Number of vertices along with ghost layers (if any) in Cartesian space. */
	float* low;			/**< Lower bound of field in continuous Cartesian space. Limit may not coincide with a grid vertex */
	float* high;			/**< Upper bound of field in continuous Cartesian space. Limit may not coincide with a grid vertex */
	int* lowInt;			/**< Lower bound of field in discrete Cartesian space (nearest grid vertex containing the lower bound). */
	int* highInt;			/**< Upper bound of field in discrete Cartesian space (nearest grid vertex containing the upper bound). */
	int* lowIntWithPad;		/**< Lower bound of field in discrete Cartesian space along with ghost layers (if any) */
	int* highIntWithPad;		/**< Upper bound of field in discrete Cartesian space along with ghost layers (if any) */

	int neighborhoodSize;   	/**< Neighborhood length for each point. Helps to determine span of ghost layers */
	int* lowPad;			/**< Span of ghost layers beyond lower end. */
	int* highPad;			/**< Span of ghost layers beyond upper end. */
	int* neighborhoodSizeArray;

public:

	/**
	 * Default constructor.
	 */
	ITL_grid(){}

	/**
	 * Constructor.
	 * @param ndim Dimensionality of the field.
	 */
	ITL_grid( int ndim )
	{
		this->dim = NULL;
		this->nDim = ndim;
	}

	/**
	 * Pure virtual function.
	 */
	virtual void setBounds(){}
	/**
	 * Pure virtual function for setting the bounds of a grid with no ghost layers.
	 */
	virtual void setBounds( float* l, float* h ){}
	/**
	 * Pure virtual function for setting the bounds of a grid with ghost layers.
	 */
	virtual void setBounds( float* l, float* h, int* lPad, int* hPad, int neighborhoodsize ){}
	virtual void setBounds( float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray ){}
	/**
	 * Pure virtual function for 3D spatial index to 1D array index conversion.
	 */
	virtual int convert3DIndex( int x, int y, int z ){ return -1; }
	/**
	 * Pure virtual function for 3D+T spatial index to 1D array index conversion.
	 */
	virtual int convert3DTimeVaryingIndex( int x, int y, int z, int t ){ return -1; }
	/**
	 * Destructor.
	 */
	virtual ~ITL_grid(){}

};

#endif
/* ITL_GRID_H_ */
