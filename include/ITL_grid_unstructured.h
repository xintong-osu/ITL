/**
 * Unstructured grid inherited from ITL_grid.
 * A class which contains information about spatial arrangement
 * and connectivity of field data points arranged without any order
 * in a Cartesian space (Not implemented).
 * Created on: Dec 9, 2010
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GRID_UNSTRUCTURED_H_
#define ITL_GRID_UNSTRUCTURED_H_

#include "ITL_grid.h"

template <class T>
class ITL_Grid_unstructured: public ITL_grid<T>
{
public:
	
	int nCell;
	
	

public:

	/**
	 * Constructor.
	 * @param ndim Dimensionality of the field.
	 */
	ITL_grid_unstructured( int ndim )
	{
		this->dim = dim;
	}
};

#endif /* ITL_GRID_UNSTRUCTURED_H_ */
