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
  // MOD-BY-LEETEN 04/09/2012-FROM:	ITL_tetrahedron() : ITL_cell(4, 4)
	ITL_tetrahedron() : ITL_cell<T>(4, 4)
  // MOD-BY-LEETEN 04/09/2012-END
	{

	}
	
};

#endif /* ITL_CELL_H_ */
