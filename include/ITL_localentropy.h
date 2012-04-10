/**
 * Local entropy field computation class.
 * Created on: Nov 19, 2010.
 * @author Abon
 * @author Teng-Yok
 * @author Cong Wang
 */

#ifndef ITL_LOCALENTROPY_H_
#define ITL_LOCALENTROPY_H_

#include "ITL_header.h"
// ADD-BY-LEETEN 07/18/2011-BEGIN
#include "ITL_histogram.h"
// ADD-BY-LEETEN 07/18/2011-END
// ADD-BY-LEETEN 04/09/2012-BEGIN
#include "ITL_grid_unstructured.h"
#include "ITL_grid_tetrahedral.h" 
#include "ITL_cell.h"
// ADD-BY-LEETEN 04/09/2012-END
#include "ITL_field_regular.h"
#include "ITL_field_unstructured.h"
#include "ITL_entropycore.h"

#include "ITL_cell.h"

template <class T>
class ITL_localentropy
{
public:

	ITL_field<T> *dataField;
	ITL_field_regular<int>* binData;	/**< A scalar field containing histogram bins corresponding to field points. */
	ITL_field_regular<float> *entropyField;

	T histogramMin;
	T histogramMax;
	bool histogramRangeSet;

	float globalEntropy;				/**< Value of computed entropy at a specified point in the field. */
	float* probarray;   				/**< Probability array used in emtropy computation. */

	ITL_histogram *histogram;			// ADD-BY-ABON 11/07/2011
	int nBin;

public:

	/**
	 * Old style Constructor.(PLAN TO REMOVE)
	 */
	ITL_localentropy( ITL_field<T> *f, ITL_histogram *hist )
	{
		this->dataField = f;
		this->binData = NULL;
		this->entropyField = NULL;
		histogramRangeSet = false;		// ADD-BY-ABON 07/19/2011
		histogram = hist;				// ADD-BY-ABON 11/07/2011
	}// End constructor

	/**
	 * Constructor.
	 */
	ITL_localentropy( ITL_field_regular<int> *bindata, ITL_histogram *hist, int nbin )
	{
		this->binData = bindata;
		histogram = hist;
		nBin = nbin;

		this->entropyField = NULL;

	}// End constructor
	
	/**
	 * Entropy computation function.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBins Number of bins used in histogram computation.
	 */
	void computeLocalEntropyOfField( bool toNormalize )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];

		// Allocate memory for entropy field (a non-padded scalar field), if not already done
		if( this->entropyField == NULL )
		{
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );


			this->entropyField = new ITL_field_regular<SCALAR>( this->dataField->getNumDim(),
																low, high );
		}

		// Compute total number of vertices in the neighborhood, including the vertext itself
		//int nNeighbors = (int)std::pow( 2.0f*this->dataField->grid->neighborhoodSize + 1.0f, this->dataField->grid->nDim );

		dataField->getNeighborhoodSize( neighborhoodSize );
		int nNeighbors = (int)( ( 2.0f*neighborhoodSize[0] + 1.0f ) *
						     ( 2.0f*neighborhoodSize[1] + 1.0f ) *
						     ( 2.0f*neighborhoodSize[2] + 1.0f ) );

		// Compute histogram if it is not already done
		//if( this->binData == NULL )
		//	this->computeHistogramBinField( nBins );

		int index1d = 0;
		int dim[4];
		entropyField->getSize( dim );

		for( int z=0; z<dim[2]; z++ )
		{
			for( int y=0; y<dim[1]; y++ )
			{
				for( int x=0; x<dim[0]; x++ )
				{
					// Compute and store the value of entropy at this point
					this->computeEntropySinglePoint( x, y, z, nNeighbors, index1d, toNormalize );

					// increment to the next element
					index1d ++;

				}// end for : x
			}// end for : y
		}// end for : z

	}// end function

	/**
	 * Entropy computation function.
	 * Creates a scalar field that contains entropy at each vertex of the unstructured grid.
	 * @param nBins Number of bins used in histogram computation.
	 */
	void computeLocalEntropyOfField_Unstructured( int nBins, bool toNormalize )
	{
		assert(dataField);
		//ITL_grid_unstructured<SCALAR>* uGrid = dynamic_cast<ITL_grid_unstructured<SCALAR>*>(dataField->grid);
		ITL_grid_unstructured<SCALAR>* uGrid = dynamic_cast<ITL_grid_unstructured<SCALAR>*>(dataField->getGrid());
		SCALAR rmin, rmax, rangeValue;

		//calculate and set histogramMin / histogramMax
		if( histogramRangeSet == false )
		{
			rmin = std::numeric_limits<double>::max();
			rmax = -std::numeric_limits<double>::max();
			for (int i = 0; i < uGrid->nCell; ++i)
			{
				ITL_cell<SCALAR>& tet = uGrid->cellList[i];
				for (int j = 0; j < tet.nVert; ++j)
				{
					double f = uGrid->vertexList[tet.v[j]].f;
					rmin = rmin < f ? rmin : f;
					rmax = rmax > f ? rmax : f;
				}
			}
			histogramMin = rmin;
			histogramMax = rmax;
			histogramRangeSet = true;
		}
		else
		{	
			rmin = histogramMin;
			rmax = histogramMax;
		}
		rangeValue = rmax - rmin;

		// Compute bin width			
		// MOD-BY-LEETEN 04/09/2012-FROM:		float binWidth = rangeValue / (float)nBin;
		float binWidth = rangeValue / (float)nBins;
		// MOD-BY-LEETEN 04/09/2012-END

		float* entropyField = new float[uGrid->getSize()];

		//compute entropy at each vertex
		for (int i = 0; i < uGrid->getSize(); ++i)
		{
			float entropy = computeEntropySinglePoint_Unstructrued(i, nBins, binWidth, toNormalize);
			entropyField[i] = entropy;
		}
	}

	/**
	 * Entropy computation function.
	 * Computes entropy at a spatial point in the field.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 * @param nNeighbors total number of neighbors in the neighborhood (including self).
	 * @param nBins Number of bins to use in histogram computation.
	 * @param entropyFieldIndex 1D index to the entopy field.
	 */
	void computeEntropySinglePoint( int x, int y, int z,
									int nNeighbors,
									int entropyFieldIndex, bool toNormalize )
	{
		#if defined( _WIN32 ) || defined( _WIN64 )
			int* binArray = new int[nNeighbors];
			int* localFreqList = new int[nBins];
		#else
			int binArray[nNeighbors];
			int localFreqList[nBin];
		#endif
		int index1d = 0;
		int binArrayIndex = 0;

		/*
		for( int k = -this->binData->grid->neighborhoodSize; k <= this->binData->grid->neighborhoodSize; k++ )
		{
			for( int j = -this->binData->grid->neighborhoodSize; j <= this->binData->grid->neighborhoodSize; j++ )
			{
				for( int i = -this->binData->grid->neighborhoodSize; i <= this->binData->grid->neighborhoodSize; i++ )
				{
		*/
		int neighborhoodSize[4];
		int lowPad[4];
		int highPad[4];
		int dimWithPad[4];
		binData->getSizeWithPad( dimWithPad );
		binData->getPadSize( lowPad, highPad );
		binData->getNeighborhoodSize( neighborhoodSize );
		//cout << dimWithPad[0] << " " << dimWithPad[1] << " " << dimWithPad[2] << " " << endl;
		//cout << lowPad[0] << " " << lowPad[1] << " " << lowPad[2] << " " << endl;
		//cout << highPad[0] << " " << highPad[1] << " " << highPad[2] << " " << endl;
		//cout << neighborhoodSize[0] << " " << neighborhoodSize[1] << " " << neighborhoodSize[2] << " " << endl;

		for( int k = -neighborhoodSize[2]; k <= neighborhoodSize[2]; k++ )
		{
			for( int j = -neighborhoodSize[1]; j <= neighborhoodSize[1]; j++ )
			{
				for( int i = -neighborhoodSize[0]; i <= neighborhoodSize[0]; i++ )
				{

					// Convert the 1D index corresponding to the point
					#ifdef DUPLICATE
					index1d = this->binData->grid->convert3DIndex( ITL_util<int>::clamp( x+i+this->binData->grid->lowPad[0], 0, this->binData->grid->dimWithPad[0]-1 ),
													ITL_util<int>::clamp( y+j+this->binData->grid->lowPad[1], 0, this->binData->grid->dimWithPad[1]-1 ),
													ITL_util<int>::clamp( z+k+this->binData->grid->lowPad[2], 0, this->binData->grid->dimWithPad[2]-1 ) );
					#endif

					#ifdef MIRROR
					index1d = this->binData->convert3DIndex( ITL_util<int>::mirror( x+i+lowPad[0], 0, dimWithPad[0]-1 ),
													ITL_util<int>::mirror( y+j+lowPad[1], 0, dimWithPad[1]-1 ),
													ITL_util<int>::mirror( z+k+lowPad[2], 0, dimWithPad[2]-1 ) );
					#endif

					// Obtain the binID corresponding to the value at this location
					//cout << i << " " << j << " " << k << " " << index1d << endl;
					binArray[binArrayIndex] = this->binData->getDataAt( index1d );

					// Move to next point in the neighborhood
					binArrayIndex ++;
				}
			}
		}

		// Compute entropy
		float entropy = ITL_entropycore::computeEntropy_HistogramBased( binArray, localFreqList, nNeighbors, nBin, toNormalize );

		// Store entropy
		this->entropyField->setDataAt( entropyFieldIndex, entropy );

		#if defined( _WIN32 ) || defined( _WIN64 )
			delete [] binArray;
			delete [] localFreqList;
		#endif

	}// end function


	/**
	 * Entropy computation function.
	 * Computes entropy at a vertex of the grid.
	 * @param vid vertex id
	 * @param nBins Number of bins to use in histogram computation.
	 * @param binWidth width of each bin
	 */
	float computeEntropySinglePoint_Unstructrued(int vid, int nBins, float binWidth, bool toNormalize)
	{
		//ITL_grid_unstructured<SCALAR>* uGrid = dynamic_cast<ITL_grid_unstructured<SCALAR>*>(dataField->grid);
		ITL_grid_unstructured<SCALAR>* uGrid = dynamic_cast<ITL_grid_unstructured<SCALAR>*>(dataField->getGrid());
		ITL_vertex<SCALAR>& vert = uGrid->vertexList[vid];

		//min/max of neighborhood box
		vec3<SCALAR> nei_min(vert.x - uGrid->radius, vert.y - uGrid->radius, vert.z - uGrid->radius);
		vec3<SCALAR> nei_max(vert.x + uGrid->radius, vert.y + uGrid->radius, vert.z + uGrid->radius);

		//get list of intersected and contained cells of vid's neighborhood box
		std::vector<int>& intersectCells = uGrid->intersectCells[vid];
		std::vector<int>& containedCells = uGrid->containCells[vid];

		//sample density
		int density = 10;

		float* hist = new float[nBins];
		memset(hist, 0, sizeof(float) * nBins);

		//each sample point's percentage in the whole cell
		float sampleFraction = 1.0/(density*density*density);

		//sample each intersected cell
		for (int i = 0; i < intersectCells.size(); ++i)
		{
			ITL_cell<SCALAR>* tet = &(uGrid->cellList[intersectCells[i]]);
			SCALAR inc = 1.0/density;
			SCALAR s, t, u;
			vec3<SCALAR> v0(uGrid->vertexList[tet->v[0]].x, uGrid->vertexList[tet->v[0]].y, uGrid->vertexList[tet->v[0]].z);
			vec3<SCALAR> v1(uGrid->vertexList[tet->v[1]].x, uGrid->vertexList[tet->v[1]].y, uGrid->vertexList[tet->v[1]].z);
			vec3<SCALAR> v2(uGrid->vertexList[tet->v[2]].x, uGrid->vertexList[tet->v[2]].y, uGrid->vertexList[tet->v[2]].z);
			vec3<SCALAR> v3(uGrid->vertexList[tet->v[3]].x, uGrid->vertexList[tet->v[3]].y, uGrid->vertexList[tet->v[3]].z);

			SCALAR f0 = uGrid->vertexList[tet->v[0]].f;
			SCALAR f1 = uGrid->vertexList[tet->v[1]].f;
			SCALAR f2 = uGrid->vertexList[tet->v[2]].f;
			SCALAR f3 = uGrid->vertexList[tet->v[3]].f;

			vec3<SCALAR> N1( v1.x - v0.x, v2.x - v0.x, v3.x - v0.x);
			vec3<SCALAR> N2( v1.y - v0.y, v2.y - v0.y, v3.y - v0.y);
			vec3<SCALAR> N3( v1.z - v0.z, v2.z - v0.z, v3.z - v0.z);
			vec3<SCALAR> vmin = nei_min - v0;
			vec3<SCALAR> vmax = nei_max - v0;

			float volume = abs(N1.dot(N2.cross(N3)) / 6.0);
			float sampleVolume = volume * sampleFraction;

			/* 
			we sample the cell in the barycentric coordinate system, s, t, u (1-s-t-u) are the coordinates for each sample,
			if a sample point is in the neighborhood box, it must satisfy
			nei_min < v1*s + v2*t + v3*u + v0*(1-s-t-u) < nei_max  in cartesian space 
			i.e.
			nei_min - v0 < (v1-v0)*s + (v2-v0)*t + (v3-v0)*u < nei_max - v0  for each x, y, z
			This defines a parallelepiped
			We then calculate the 8 corners of the parallelepiped, so to get the range of s, t, u, hopefully they will be much 
			smaller than [0,1]. Thus we can do less sampling.

			Corner of parallelepiped is given by the intersection of 3 planes
		    N1 . p = d1
		    N2 . p = d2
		    N3 . p = d3
		    In the above and what follows, "." signifies the dot product and "*" is the cross product. 
		    The intersection point P is given by:
		                 d1 ( N2 * N3 ) + d2 ( N3 * N1 ) + d3 ( N1 * N2 )
		    P = 	-------------------------------------------------------------------------
		                               N1 . ( N2 * N3 ) 
			*/

			vec3<SCALAR> C23 = N2.cross(N3);
			vec3<SCALAR> C31 = N3.cross(N1);
			vec3<SCALAR> C12 = N1.cross(N2);
			SCALAR D1C23 = N1.dot(C23);
			SCALAR stumin[3], stumax[3];
			if (D1C23)
			{
				//get 8 corners
				vec3<SCALAR> p[8];
				/*
				p[0] = ( C23 * vmin.x + C31 * vmin.y + C12 * vmin.z) / D1C23;
				p[1] = ( C23 * vmin.x + C31 * vmin.y + C12 * vmax.z) / D1C23;
				p[2] = ( C23 * vmin.x + C31 * vmax.y + C12 * vmin.z) / D1C23;
				p[3] = ( C23 * vmin.x + C31 * vmax.y + C12 * vmax.z) / D1C23;
				p[4] = ( C23 * vmax.x + C31 * vmin.y + C12 * vmin.z) / D1C23;
				p[5] = ( C23 * vmax.x + C31 * vmin.y + C12 * vmax.z) / D1C23;
				p[6] = ( C23 * vmax.x + C31 * vmax.y + C12 * vmin.z) / D1C23;
				p[7] = ( C23 * vmax.x + C31 * vmax.y + C12 * vmax.z) / D1C23;
				*/

				//get s,t,u range
				stumin[0] = stumin[1] = stumin[2] = 1; 
				stumax[0] = stumax[1] = stumax[2] = 0; 
				for (int i = 0; i < 8; ++i)
				{
					stumin[0] = stumin[0] < p[i].x ? stumin[0] : p[i].x;
					stumin[1] = stumin[1] < p[i].y ? stumin[1] : p[i].y;
					stumin[2] = stumin[2] < p[i].z ? stumin[2] : p[i].z;
					stumax[0] = stumax[0] > p[i].x ? stumax[0] : p[i].x;
					stumax[1] = stumax[1] > p[i].y ? stumax[1] : p[i].y;
					stumax[2] = stumax[2] > p[i].z ? stumax[2] : p[i].z;
				}
				stumin[0] = stumin[0] > 0 ? stumin[0] : 0;
				stumin[1] = stumin[1] > 0 ? stumin[1] : 0;
				stumin[2] = stumin[2] > 0 ? stumin[2] : 0;
				stumax[0] = stumax[0] < 1 ? stumax[0] : 1;
				stumax[1] = stumax[1] < 1 ? stumax[1] : 1;
				stumax[2] = stumax[2] < 1 ? stumax[2] : 1;
			}
			else
			{
				stumin[0] = stumin[1] = stumin[2] = 0; 
				stumax[0] = stumax[1] = stumax[2] = 1;
			}

			//new starting point for s, t, u 
			SCALAR s0 = (int(stumin[0] / inc))*inc;
			SCALAR t0 = (int(stumin[1] / inc))*inc;
			SCALAR u0 = (int(stumin[2] / inc))*inc;

			//doing incremental calculation to avoid expensive multiplication
			vec3<SCALAR> ds((v1-v0)*inc), dt((v2-v0)*inc), du((v3-v0)*inc);
			vec3<SCALAR> cs, ct, cu;
			vec3<SCALAR> cs0((v1-v0)*s0 - ds), ct0((v2-v0)*t0 - dt), cu0((v3-v0)*u0 -du);

			SCALAR dfs((f1-f0)*inc), dft((f2-f0)*inc), dfu((f3-f0)*inc);
			SCALAR cfs, cft, cfu;
			SCALAR cfs0((f1-f0)*s0 - dfs), cft0((f2-f0)*t0 - dft), cfu0((f3-f0)*u0 -dfu);

			cs = cs0;
			cfs = cfs0;
			for (s = s0; s <= stumax[0]; s += inc)
			{
				cs = cs + ds;
				ct = cs + ct0;
				cfs = cfs + dfs;
				cft = cfs + cft0;
				for (t = t0; t <= stumax[1] && s + t <=1.00001; t += inc)
				{
					ct = ct + dt;
					cu = ct + cu0;
					cft = cft + dft;
					cfu = cft + cfu0;
					for (u = u0; u <= stumax[2] && s + t + u <= 1.00001; u += inc)
					{
						cu = cu + du;
						cfu = cfu + dfu;
						vec3<SCALAR> v = cu + v0;

						//if sample points is in box, add contribution to histogram
						if (v.x >= nei_min.x && v.y >= nei_min.y 
							&& v.z >= nei_min.z && v.x <= nei_max.x && v.y <= nei_max.y && v.z <= nei_max.z)
						{
							SCALAR f = cfu + f0;
							int binId = (int)floor( ( f - histogramMin ) / binWidth  );
							binId = ITL_util<int>::clamp( binId, 0, nBins-1 );
							hist[binId] += sampleVolume;
						}
					}
				}
			}
		}

		//for contained cells, do contour spectrum
		double func[21]= {0};
		double min=0,max=0;
		int posmin=0,posmax=0;
		double tt = 0;
		ITL_grid_tetrahedral<SCALAR>* tetGrid = dynamic_cast<ITL_grid_tetrahedral<SCALAR>*>(uGrid);
		for (int i = 0; i < containedCells.size(); ++i)
		{
			tetGrid->local_func(i,func);

			//min and max defines range of cell i
			min=func[0];
			max=func[5];
			posmin = floor((min - histogramMin) / binWidth);
			for(int j = posmin; histogramMin + j * binWidth < max; j++)
			{
				tt = tetGrid->contribution(func,histogramMin + j * binWidth, histogramMin + (j + 1) * binWidth);
				if(tt > 0 || tt < 0)
					hist[j] += tt;
			}
		}

		float entropy = ITL_entropycore::computeEntropy_HistogramBased( hist, nBins, toNormalize );
		return entropy;
	}

	// ADD-BY-ABON 07/19/2011-BEGIN	
	/**
	  * Historgam range set function
	  * Sets the range over which historgam computation is to be done.
	  * @param minR lower limit of the range.
	  * @param maxR upper limit of the range.
	  */
	void setHistogramRange( T minR, T maxR )
	{
		histogramMin = minR;
		histogramMax = maxR;
		histogramRangeSet = true;
	}
	// ADD-BY-ABON 07/19/2011-END	


	/**
	 * Entropy accessor function.
	 * @return computed entropy.
	 */
	float getEntropy();

	/**
	 * Entropy field accessor function.
	 * Returns computed entropy field.
	 * @return pointer to entropy field.
	 */
	ITL_field_regular<float>* getEntropyField()
	{
		return this->entropyField;
	}// end function
};

#endif
/* ITL_LOCALENTROPY_H_ */
