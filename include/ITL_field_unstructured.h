/**
 * Unstructured field class inherited from ITL_field.
 * Container for irregularly spaced static/time-varying scalar/vector data.
 * Created on: Dec 09, 2010
 * @author Abon
 * @author Teng-Yok
 * @author Cong Wang
 */

#ifndef ITL_FIELD_UNSTRUCTURED_H_
#define ITL_FIELD_UNSTRUCTURED_H_

#include "ITL_grid_tetrahedral.h"

template <class T>
class ITL_field_unstructured: public ITL_field<T>
{
  using ITL_field<T>::grid; // ADD-BY-LEETEN 04/09/2012
public:

	ITL_field_unstructured(ITL_grid_unstructured<T>* ugrid)
	{
		grid = ugrid;
		this->datastore = new ITL_datastore<T>(grid->nVertices);
	}
};

#endif /* ITL_FIELD_UNSTRUCTURED_H_ */
