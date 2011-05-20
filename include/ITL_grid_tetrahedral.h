/**
 * Tetrahedral unstructured grid inherited from ITL_grid_unstructured.
 * A class which contains information about spatial arrangement
 * and connectivity of field data points arranged as a series of tetrahedra
 * in a Cartesian space. We currently assume that the grid is 3-dimensional.
 * Created on: Aptil 28, 2011
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GRID_TETRAHEDRAL_H_
#define ITL_GRID_TETRAHEDRAL_H_

#include "ITL_grid.h"
#include "ITL_grid_unstructured.h"
#include "ITL_tetrahedron.h"

template <class T>
class ITL_grid_tetrahedral: public ITL_grid_unstructured<T>
{
public:
	
	int nCell;
	ITL_tetrahedron *listOfCells;

public:

	/**
	 * Constructor.
	 * @param ndim Dimensionality of the field.
	 */
	ITL_grid_tetrahedral()
	{
		this->dim = dim;
		this->listOfCells = NULL;		
	}
};

#endif /* ITL_GRID_TETRAHEDRAL_H_ */
