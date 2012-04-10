/**
 * Tetrahedral unstructured grid inherited from ITL_grid_unstructured.
 * A class which contains information about spatial arrangement
 * and connectivity of field data points arranged as a series of tetrahedra
 * in a Cartesian space. We currently assume that the grid is 3-dimensional.
 * Created on: Aptil 28, 2011
 * @author Abon
 * @author Teng-Yok
 * @author Cong Wang
 * @author Tiantian Xu
 */

#ifndef ITL_GRID_TETRAHEDRAL_H_
#define ITL_GRID_TETRAHEDRAL_H_


#include "ITL_grid_unstructured.h"
#include "ITL_tetrahedron.h"

#define EFFICIENT

template <class T>
class ITL_grid_tetrahedral: public ITL_grid_unstructured<T>
{
  // MOD-BY-LEETEN 04/09/2012-FROM: private:
public:
  using ITL_grid_unstructured<T>::nDim;
  using ITL_grid_unstructured<T>::radius; 
  using ITL_grid_unstructured<T>::nVertices;
  using ITL_grid_unstructured<T>::nCell;
  using ITL_grid_unstructured<T>::vertexList; //list of vertices in the grid
  using ITL_grid_unstructured<T>::cellList; //list of cells in the grid
  using ITL_grid_unstructured<T>::intersectCells; //vector of intersecting cells of each vertex's neighborhood box in the grid
  using ITL_grid_unstructured<T>::containCells; //vector of contained cells of each vertex's neighborhood box in the grid
  // MOD-BY-LEETEN 04/09/2012-END
	//check if bbox intersects local entropy box
	inline bool boxIntersec(const ITL_vertex<T>& v, const bbox<T>& bb)
	{
		bool isIntersec = true;
		if ((bb.min.x >= v.x + radius) 
			|| (bb.min.y >= v.y + radius) || (bb.min.z >= v.z + radius) 
			|| (bb.max.x <= v.x - radius) || (bb.max.y <= v.y - radius)
			|| (bb.max.z <= v.z - radius))
		{
			isIntersec = false;
		}
		return isIntersec;
	}
	
	//check if bbox is contained in local entropy box
	inline bool boxContained(const ITL_vertex<T>& v, const bbox<T>& bb)
	{
		bool isContained = false;
		if ((bb.max.x <= v.x + radius) 
			&& (bb.max.y <= v.y + radius) && (bb.max.z <= v.z + radius) 
			&& (bb.min.x >= v.x - radius) && (bb.min.y >= v.y - radius)
			&& (bb.min.z >= v.z - radius))
		{
			isContained = true;
		}
		return isContained;
	}

	//seed propagation algorithm for each vertex to find its neighborhood box's intersected cells and contained cells
	void seedPropagation(const int verIndex, int* tetVisitedHash, vector<int>& intersectTets, vector<int>& containTets)
	{
		const ITL_vertex<T>& v = vertexList[verIndex];
		deque<int> possibleIntersecTets;

		possibleIntersecTets.insert(possibleIntersecTets.end(), v.adjacentCells->begin(), v.adjacentCells->end());
		for (vector<int>::iterator it = v.adjacentCells->begin(); it < v.adjacentCells->end(); ++it)
		{
			tetVisitedHash[*it] = verIndex;
		}
		while(!possibleIntersecTets.empty())
		{
			int checkIndex = possibleIntersecTets.front();
			const bbox<T>& bb = cellList[checkIndex].getBBox();        
			if (boxIntersec(v, bb))
			{
				if (!boxContained(v, bb))
				{
					intersectTets.push_back(checkIndex);
				}
				else
				{
					containTets.push_back(checkIndex);
				}

				for (int i = 0; i < 4; ++i)
				{
					int index = cellList[checkIndex].neighborCell[i];
					if (index != -1 && tetVisitedHash[index] != verIndex)
					{
						possibleIntersecTets.push_back(index);
						tetVisitedHash[index] = verIndex;
					}                
				}            
			}
			possibleIntersecTets.pop_front();
		}
	}

public:

	/**
	 * Constructor.
	 * @param nV # of vertices
	 * @param nC # of cells
	 */
	ITL_grid_tetrahedral(int nV, int nC)
	{
		this->nDim = 3;
		this->nVertices = nV;
		this->nCell = nC;
		vertexList = new ITL_vertex<T>[nVertices];
		cellList = new ITL_tetrahedron<T>[nCell];	
		this->intersectCells = new vector<int>[this->nVertices];
		this->containCells = new vector<int>[this->nVertices];
	}

	ITL_tetrahedron<T>& getTet(int index)
	{
		return cellList[index];
	}

	//building tetrahedron adjacent information
	void buildTetAdjacentInfor()
	{
		//building vertex adjacent information
		for (int i = 0; i < nCell; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				this->vertexList[this->cellList[i].v[j]].adjacentCells->push_back(i);
			}
		}

		//building tetrahedron adjacent information
		for (int k = 0; k < nCell; ++k)
		{
			for (int i = 0; i < 4; ++i)
			{
				ITL_cell<SCALAR>* tet = &this->cellList[k];
				std::set<int> hash1, hash2;
				int j = i;
				const ITL_vertex<SCALAR>& a = this->vertexList[tet->v[++j % 4]];
				const ITL_vertex<SCALAR>& b = this->vertexList[tet->v[++j % 4]];
				const ITL_vertex<SCALAR>& c = this->vertexList[tet->v[++j % 4]];
				for (vector<int>::iterator it = a.adjacentCells->begin(); it < a.adjacentCells->end(); ++it)
				{
					hash1.insert(*it);
				}
				for (vector<int>::iterator it = b.adjacentCells->begin(); it < b.adjacentCells->end(); ++it)
				{
					if (hash1.find(*it) != hash1.end())
					{
						hash2.insert(*it);
					}
				}
				for (vector<int>::iterator it = c.adjacentCells->begin(); it < c.adjacentCells->end(); ++it)
				{
					if (hash2.find(*it) != hash2.end() && *it != tet->index)
					{
						tet->neighborCell[i] = *it;
						break;
					}
				}
			}
		}
	}

	//build intersected and contained cells information, possibly with external vectors
	void getBoxIntersecTets(vector<int>*& intersectTets, vector<int>*& containTets)
	{
		cout << "start BoxIntersecTets test: " << endl;
		int* tetVisitedHash = new int [nCell];
		memset(tetVisitedHash, -1, sizeof(int)*nCell);
		for (int i = 0; i < nVertices; ++i)
		{
			seedPropagation(i,tetVisitedHash, intersectTets[i], containTets[i]);
		}
		delete[] tetVisitedHash;
		cout << "finished." << endl;
	}

	//build intersected and contained cells information, with internal vectors (intersectCells / containCells)
	void buildBoxIntersecTets()
	{
		getBoxIntersecTets(this->intersectCells, this->containCells);
	}

	void sortVerts(double arr[],int n, int *pos)
	{
		for (int i = 1; i < n; ++i)
		{
			for (int j = i - 1; arr[j] > arr[j+1] && j >= 0; --j)
			{
				double temp = arr[j];
				arr[j] = arr[j+1];
				arr[j+1] = temp;	

				int t=pos[j];
				pos[j]=pos[j+1];
				pos[j+1]=t;
			}			
		}		
	}

	void linterp(double s1,double s2,double s,double p1[],double p2[],double *p)
	{
		double ratio=(s-s1)/(s2-s1);

		for(int i=0; i<3; i++)
			p[i]=p1[i]+ratio*(p2[i]-p1[i]);
	}

	//calculate area of the triangle defined by p1,p2 and p3
	double triangle_area(double *p1,double *p2, double *p3)
	{
		double AB=(p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]);
		double BC=(p1[0]-p3[0])*(p1[0]-p3[0])+(p1[1]-p3[1])*(p1[1]-p3[1])+(p1[2]-p3[2])*(p1[2]-p3[2]);
		double CA=(p3[0]-p2[0])*(p3[0]-p2[0])+(p3[1]-p2[1])*(p3[1]-p2[1])+(p3[2]-p2[2])*(p3[2]-p2[2]);

		AB=sqrt(AB);
		BC=sqrt(BC);
		CA=sqrt(CA);

		double p=(AB+BC+CA)/2;

		double area=p*(p-AB)*(p-BC)*(p-CA);
		area=sqrt(area);

		return area;
	}

	double cmptl(int cid, int tgtid)
	{
		double area=0;

		//retrieve from cell list the vertices in current triangle
		int v1,v2,v3,v4;
		v1=cellList[cid].v[0];
		v2=cellList[cid].v[1];
		v3=cellList[cid].v[2];
		v4=cellList[cid].v[3];

		double list[4]= {vertexList[v1].f,vertexList[v2].f,vertexList[v3].f,vertexList[v4].f};
		int pos[4]= {v1,v2,v3,v4};
		sortVerts(list,4,pos);//pos stores the index of sorted elements, list[pos[0]]<list[pos[1]]<...

		double p1[3]= {vertexList[pos[0]].x,vertexList[pos[0]].y,vertexList[pos[0]].z};
		double p2[3]= {vertexList[pos[1]].x,vertexList[pos[1]].y,vertexList[pos[1]].z};
		double p3[3]= {vertexList[pos[2]].x,vertexList[pos[2]].y,vertexList[pos[2]].z};
		double p4[3]= {vertexList[pos[3]].x,vertexList[pos[3]].y,vertexList[pos[3]].z};

		double pi11[3],pi21[3];
		//Do linear interpolation and compute intersection area at scalar value=list[1]
		if(tgtid==1)
		{
			linterp(list[0],list[2],list[1],p1,p3,pi11);
			linterp(list[0],list[3],list[1],p1,p4,pi21);
			area=triangle_area(pi11,pi21,p2);
			return area;

		}
		//Do linear interpolation and compute intersection area at scalar value=list[2]
		else if(tgtid==2)
		{
			linterp(list[1],list[3],list[2],p2,p4,pi11);
			linterp(list[0],list[3],list[2],p1,p4,pi21);
			area=triangle_area(pi11,pi21,p3);
			return area;
		}

		return area;
	}


	//The 'efbspline' is an efficient way to do spline computation, a result of mathematical equivalence
	//An observation was made that the final spline functions, on the range of t[0]-t[5] consists of three, not five different quadratic functions
	//The first part is t[0]-t[2], second part is t[2]-t[3], third part is t[3]-t[5]
	//such obeservation is consistent with the tetrahedra's geometric property (in terms of the volume of the tetra)
	//Here we save some time by not computing those basis functions that we don't need
	void efbspline(double t[],vec3<SCALAR> cp[],double *f1,double *f2,double *f3)
	{
		double b11[3],b13[3],b22[3],b31[3],b33[3];
		double a=0,b=0,c=0;
		double m=0,n=0;

		//for basis b-spline function on t0-t3,only compute first and third part
		//part 1:t0<=t<t1
		if(t[2]!=t[0])
			m=1/(t[2]-t[0]);
		if(t[1]!=t[0])
			n=1/(t[1]-t[0]);

		a=m*n;
		b=-2*t[0]*m*n;
		c=t[0]*t[0]*m*n;

		b11[0]=a;
		b11[1]=b;
		b11[2]=c;

		//part 3: t2<=t<t3
		m=0,n=0;
		if(t[3]!=t[1])
			m=1/(t[3]-t[1]);
		if(t[3]!=t[2])
			n=1/(t[3]-t[2]);

		a=m*n;
		b=-2*m*n*t[3];
		c=m*n*t[3]*t[3];

		b13[0]=a;
		b13[1]=b;
		b13[2]=c;

		//for basis b-spline function on t1-t4, compute the second part
		//t2<=t<=t3
		double u=0,v=0,w=0;
		if(t[3]!=t[1])
			u=1/(t[3]-t[1]);
		if(t[3]!=t[2])
			v=1/(t[3]-t[2]);
		if(t[4]!=t[2])
			w=1/(t[4]-t[2]);

		a=-(u*v+w*v);
		b=u*v*(t[1]+t[3])+w*v*(t[2]+t[4]);
		c=-u*v*t[1]*t[3]-v*t[2]-w*v*t[2]*t[2];

		b22[0]=a;
		b22[1]=b;
		b22[2]=c;

		//for basis b-spline function on t2-t5, compute the first and third part
		//part 1:t2<=t<t3
		if(t[4]!=t[2])
			m=1/(t[4]-t[2]);
		if(t[3]!=t[2])
			n=1/(t[3]-t[2]);

		a=m*n;
		b=-2*t[2]*m*n;
		c=t[2]*t[2]*m*n;

		b31[0]=a;
		b31[1]=b;
		b31[2]=c;

		//part 3: t4<=t<t5
		m=0,n=0;
		if(t[5]!=t[3])
			m=1/(t[5]-t[3]);
		if(t[5]!=t[4])
			n=1/(t[5]-t[4]);

		a=m*n;
		b=-2*m*n*t[5];
		c=m*n*t[5]*t[5];

		b33[0]=a;
		b33[1]=b;
		b33[2]=c;

		//compute b-spline function
		//(t0,t2): y=b11*cp[0].y
		f1[0]=b11[0]*cp[0].y;
		f1[1]=b11[1]*cp[0].y;
		f1[2]=b11[2]*cp[0].y;

		//(t2,t3): y=b13*cp[0].y+b22*cp[1].y+b31*cp[2].y
		f2[0]=b13[0]*cp[0].y+b22[0]*cp[1].y+b31[0]*cp[2].y;
		f2[1]=b13[1]*cp[0].y+b22[1]*cp[1].y+b31[1]*cp[2].y;
		f2[2]=b13[2]*cp[0].y+b22[2]*cp[1].y+b31[2]*cp[2].y;

		//(t3,t5): y=b33*cp[2].y
		f3[0]=b33[0]*cp[2].y;
		f3[1]=b33[1]*cp[2].y;
		f3[2]=b33[2]*cp[2].y;
	}

	void local_func(int cid,double *func)
	{
		//Compute the intersection area with the contour and the tetrahedra at scalar=second smallest scalar in cell 'cid'
		double areav2=cmptl(cid,1);
		//Compute the intersection area with the contour and the tetrahedra at scalar=second largest scalar in cell 'cid'
		double areav3=cmptl(cid,2);

		int v1,v2,v3,v4;
		v1=cellList[cid].v[0];
		v2=cellList[cid].v[1];
		v3=cellList[cid].v[2];
		v4=cellList[cid].v[3];

		double list[4]= {vertexList[v1].f,vertexList[v2].f,vertexList[v3].f,vertexList[v4].f};
		int pos[4]= {v1,v2,v3,v4};
		sortVerts(list,4,pos);//pos stores the index of sorted elements, list[pos[0]]<list[pos[1]]<...

		//6 knots, 3 control points, degree 2
		double t[6]= {0};

		t[0]=list[0];
		t[1]=(list[0]+list[1])/2;
		t[2]=list[1];
		t[3]=list[2];
		t[4]=(list[2]+list[3])/2;
		t[5]=list[3];

		double x=areav2/(list[1]-list[0]);
		double y=areav3/(list[2]-list[3]);

		double px=(list[1]+list[2])/2.0;
		double py=areav2*(list[2]-list[0])/(list[1]-list[0]);

		double value3=py*(list[3]-list[2])/(list[3]-list[1]);

		vec3<SCALAR> cp[3];
		cp[0].x=t[2];
		cp[0].y=areav2/2.0;
		cp[0].z=0;

		cp[1].x=px;
		cp[1].y=py;
		cp[1].z=0;

		cp[2].x=t[3];
		cp[2].y=value3/2.0;
		cp[2].z=0;

		double f1[3],f2[3],f3[3];

		efbspline(t,cp,f1,f2,f3);

		int i=0;

		for(i=0; i<6; i++)
			func[i]=t[i];

		for(i=6; i<9; i++)
			func[i]=f1[i-6];

		for(i=9; i<12; i++)
			func[i]=f2[i-9];

		for(i=12; i<15; i++)
			func[i]=f3[i-12];
	}

	//Integration of a quadratic function 'f' on [min, max]
	double integr(double f[],double min, double max)
	{
		if(min==max||min>max) return 0;

		double result=0;

		double p3=0,p2=0,p1=0;

		p3=f[0]/3; 
		p2=f[1]/2;
		p1=f[2];

		result=p3*(pow(max,3)-pow(min,3))+p2*(pow(max,2)-pow(min,2))+p1*(max-min);

		return result;
	}

	//Compute integration of a piecewise quadratic function 'func' on the range [min,max]
	double contribution(double func[],double min, double max)
	{
		int i=0,p=0;
		double v[4],f1[3],f2[3],f3[3];
		double result=0;
		v[0]=func[0];
		v[1]=func[2];
		v[2]=func[3];
		v[3]=func[5];

		for(i=0;i<3;i++) f1[i]=func[i+6];
		for(i=0;i<3;i++) f2[i]=func[i+9];
		for(i=0;i<3;i++) f3[i]=func[i+12];

		if(min<=v[0])
		{
			if(max<=v[0]) return 0;
			else if(max<=v[1])
				result+=integr(f1,v[0],max);
			else if(max<=v[2])
			{
				result+=integr(f1,v[0],v[1]);
				result+=integr(f2,v[1],max);
			}else if(max<=v[3])
			{
				result+=integr(f1,v[0],v[1]);
				result+=integr(f2,v[1],v[2]);
				result+=integr(f3,v[2],max);

			}else if(max>v[3])
			{
				result+=integr(f1,v[0],v[1]);
				result+=integr(f2,v[1],v[2]);
				result+=integr(f3,v[2],v[3]);

			}else return 0;
			return result;
		}
		else if(min<=v[1])
		{
			if(max<=v[1])
				result+=integr(f1,min,max);
			else if(max<=v[2])
			{
				result+=integr(f1,min,v[1]);
				result+=integr(f2,v[1],max);
			}else if(max<=v[3])
			{
				result+=integr(f1,min,v[1]);
				result+=integr(f2,v[1],v[2]);
				result+=integr(f3,v[2],max);

			}else if(max>v[3])
			{
				result+=integr(f1,min,v[1]);
				result+=integr(f2,v[1],v[2]);
				result+=integr(f3,v[2],v[3]);
			}else return 0;

			return result;
		}
		else if(min<=v[2])
		{
			if(max<=v[2])
			{
				result+=integr(f2,min,max);
			}else if(max<=v[3])
			{
				result+=integr(f2,min,v[2]);
				result+=integr(f3,v[2],max);

			}else if(max>v[3])
			{
				result+=integr(f2,min,v[2]);
				result+=integr(f3,v[2],v[3]);

			}else return 0;
			return result;
		}
		else if(min<=v[3])
		{
			if(max<=v[3])
			{
				result+=integr(f3,min,max);

			}else  if(max>v[3])
			{
				result+=integr(f3,min,v[3]);

			}else return 0;
			return result;
		}else return 0;
	}


	double glb_entropy(int hisize)
	{
		int i=0,j=0;
		//extract scalar list
		double *scalar_list=new double[nVertices];
		for(i=0; i<nVertices; i++)
			scalar_list[i]=vertexList[i].f;

		double rangemin=scalar_list[0],rangemax=scalar_list[0];
		for(i=1; i<nVertices; i++)
		{
			if(scalar_list[i]>rangemax) rangemax=scalar_list[i];
			if(scalar_list[i]<rangemin) rangemin=scalar_list[i];
		}

		double list[1001]= {0};
		double delta=(rangemax-rangemin)/hisize;//the span of every bin
		for(i=0; i<hisize+1; i++)
		{
			list[i]=rangemin+i*delta;//hisize bins, bin i starts with scalar value list[i]
		}

		//process cell by cell: calculate the contribution of cell i to bin j
		double min=0,max=0;
		int posmin=0,posmax=0;

		double histogram[1001]= {0};
		double func[21]= {0};
		double tt=0;

		for(i=0; i<nCell; i++)
		{
			local_func(i,func);

			//min and max defines range of cell i
			min=func[0];
			max=func[5];

			posmin=floor((min-rangemin)/delta);
			for(j=posmin; list[j]<max; j++)
			{
				tt=contribution(func,list[j],list[j+1]);
				if(tt>=0||tt<=0)
					histogram[j]+=tt;
			}
		}

		//Compute the probabilities
		double total=0;
		for(i=0; i<hisize; i++)
		{
			total+=histogram[i];
		}

		double prb[1000]= {0};
		for(i=0; i<hisize; i++)
			prb[i]=histogram[i]/total;

		double entropy=0;
		//Compute the entropy	
		for(i=0; i<hisize; i++)
		{
			if(prb[i]>0)
				entropy-=prb[i]*(log(prb[i])/log(2.0));
		}
		return entropy;
	}
};

#endif /* ITL_GRID_TETRAHEDRAL_H_ */
