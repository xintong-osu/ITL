/**
 * Tetrahedron class inherited from ITL_cell.
 * Encapsulates a tetrahedral unit of an unstructured grid.
 * Created on: Apr 28, 2011
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_TETRAHEDRON_H_
#define ITL_TETRAHEDRON_H_

#include "ITL_cell.h"

template <class T>
class ITL_tetrahedron: public ITL_cell<T>
{
public:
	ITL_tetrahedron() : ITL_cell(4, 4)
	{

	}
	
};

#endif /* ITL_CELL_H_ */
