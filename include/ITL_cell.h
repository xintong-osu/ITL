/**
 * Cell class representing a unit of an unstructured grid.
 * Created on: Apr 28, 2011
 * @author Abon
 * @author Teng-Yok
 * @author Cong Wang
 */

#ifndef ITL_CELL_H_
#define ITL_CELL_H_

#include "ITL_vectormatrix.h"
#include <limits>

//vertex class used in unstructured grid 
template <class T>
struct ITL_vertex 
{
	float x;
	float y;
	float z;
	T f;
	vector<int>* adjacentCells;  //only needed when build tet adjacency infor, get rid of later

	ITL_vertex()
	{
		x = y = z = f = 0;
		adjacentCells = new vector<int>;
	}

	ITL_vertex(const ITL_vertex& right)
	{
		x = right.x;
		y = right.y;
		z = right.z;
		f = right.f;
		adjacentCells = new vector<int>;
		*adjacentCells = *right.adjacentCells;
	}

	~ITL_vertex()
	{
		delete adjacentCells;
	}
};

template <class T>
struct vec3
{
	T x;
	T y;
	T z;

	vec3():x(0), y(0), z(0)
	{
	}

	//vec3(vertex& right):x(right.x), y(right.y), z(right.z)
	//{
	//}

	vec3(T rx, T ry, T rz):x(rx),y(ry),z(rz)
	{
	}

#if 0 // MOD-BY-LEETEN 04/09/2012-FROM:
	vec3 cross(const vec3& r)
	{
		return vec3(this->y * r.z - this->z * r.y, this->z * r.x - this->x * r.z, this->x * r.y - this->y * r.x);
	}

	T dot(const vec3& r)
	{
		return this->x * r.x + this->y * r.y + this->z * r.z;
	}

	vec3<T> operator- (vec3<T>& right)
	{
		vec3<T> v(x-right.x, y-right.y, z-right.z);
		return v;
	}

	vec3 operator+ (vec3& right)
	{
		vec3 v(x+right.x, y+right.y, z+right.z);
		return v;
	}

	vec3 operator* (T f)
	{
		vec3 v(x*f, y*f, z*f);
		return v;
	}

	vec3 operator/ (T f)
	{
		vec3 v(x/f, y/f, z/f);
		return v;
	}
#else // MOD-BY-LEETEN 04/09/2012-TO:
	const vec3 cross(const vec3& r) const
	{
		return vec3(this->y * r.z - this->z * r.y, this->z * r.x - this->x * r.z, this->x * r.y - this->y * r.x);
	}

	T dot(const vec3& r) const
	{
		return this->x * r.x + this->y * r.y + this->z * r.z;
	}

	const vec3<T>& operator- (const vec3<T>& right) const
	{
		vec3<T> v(x-right.x, y-right.y, z-right.z);
		return v;
	}

	const vec3<T>& operator+ (const vec3& right) const
	{
		vec3 v(x+right.x, y+right.y, z+right.z);
		return v;
	}

	const vec3<T>& operator* (T f) const
	{
		vec3 v(x*f, y*f, z*f);
		return v;
	}

	const vec3<T>& operator/ (T f) const
	{
		vec3 v(x/f, y/f, z/f);
		return v;
	}
#endif // MOD-BY-LEETEN 04/09/2012-END

	vec3& operator*= (T f)
	{
		x = x*f;
		y = y*f;
		z = z*f;
		return *this;
	}

	vec3(const vec3& right):x(right.x),y(right.y),z(right.z)
	{
	}

};

//bounding box for cells
template <class T>
struct bbox
{
	vec3<T> min;
	vec3<T> max;
	bbox()
	{
		min.x = min.y = min.z = numeric_limits<T>::max();
		max.x = max.y = max.z = -numeric_limits<T>::max();
	}
};

template <class T>
class ITL_cell
{
public:
	int index;
	int nVert;
	int* v;
	int nFace;
	int* neighborCell;
	bbox<T> bb;

	ITL_cell(int nv, int nf)
	{
		nVert = nv;
		nFace = nf;

		index = -1;
		v = new int[nVert];
		for (int i = 0; i < nVert; ++i)
		{
			v[i] = -1;
		}
		neighborCell = new int[nFace];
		for (int i = 0; i < nFace; ++i)
		{
			neighborCell[i] = -1;
		}
	}

	~ITL_cell()
	{
		delete[] v;
		delete[] neighborCell;
	}

	bbox<T>& getBBox()
	{
		return bb;
	}

	void setBBox(T minx, T miny, T minz, T maxx, T maxy, T maxz)
	{
		bb.min.x = minx;
		bb.min.y = miny;
		bb.min.z = minz;
		bb.max.x = maxx;
		bb.max.y = maxy;
		bb.max.z = maxz;
	}
};

#endif /* ITL_CELL_H_ */
