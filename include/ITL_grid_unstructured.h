/**
 * Unstructured grid inherited from ITL_grid.
 * A class which contains information about spatial arrangement
 * and connectivity of field data points arranged without any order
 * in a Cartesian space (Not implemented).
 * Created on: Dec 9, 2010
 * @author Abon
 * @author Teng-Yok
 * @author Cong Wang
 */

#ifndef ITL_GRID_UNSTRUCTURED_H_
#define ITL_GRID_UNSTRUCTURED_H_

#include "ITL_grid.h"
#include "ITL_cell.h"
#include <vector>

template <class T>
class ITL_grid_unstructured: public ITL_grid<T>
{
public:
	int nCell;
	ITL_vertex<T>* vertexList; //list of vertices in the grid
	ITL_cell<T>* cellList; //list of cells in the grid
	
	std::vector<int>* intersectCells; //vector of intersecting cells of each vertex's neighborhood box in the grid
	std::vector<int>* containCells; //vector of contained cells of each vertex's neighborhood box in the grid

	T radius; //size of neighborhood box

public:

	//build bounding box information for each cell
	void buildBoundingBoxInfor()
	{
		for (int i = 0; i < nCell; ++i)
		{
			double minx, miny, minz, maxx, maxy, maxz;
			minx = miny = minz = std::numeric_limits<double>::max();
			maxx = maxy = maxz = -std::numeric_limits<double>::max();
			for (int j = 0; j < 4; ++j)
			{
				maxx = this->vertexList[this->cellList[i].v[j]].x > maxx ? this->vertexList[this->cellList[i].v[j]].x : maxx;
				maxy = this->vertexList[this->cellList[i].v[j]].y > maxy ? this->vertexList[this->cellList[i].v[j]].y : maxy;
				maxz = this->vertexList[this->cellList[i].v[j]].z > maxz ? this->vertexList[this->cellList[i].v[j]].z : maxz;
				minx = this->vertexList[this->cellList[i].v[j]].x < minx ? this->vertexList[this->cellList[i].v[j]].x : minx;
				miny = this->vertexList[this->cellList[i].v[j]].y < miny ? this->vertexList[this->cellList[i].v[j]].y : miny;
				minz = this->vertexList[this->cellList[i].v[j]].z < minz ? this->vertexList[this->cellList[i].v[j]].z : minz;
			}
			this->cellList[i].setBBox(minx,  miny,  minz,  maxx,  maxy,  maxz);
		}
	}

};

#endif /* ITL_GRID_UNSTRUCTURED_H_ */
