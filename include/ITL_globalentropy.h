/**
 * Gloabl entropy computation class.
 * Created on: April 02, 2011.
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GLOBALENTROPY_H_
#define ITL_GLOBALENTROPY_H_

#include "ITL_header.h"
#include "ITL_field_regular.h"
#include "ITL_histogram.h"
#include "ITL_histogrammapper.h"
#include "ITL_entropycore.h"

template <class T>
class ITL_globalentropy
{
public:

	ITL_field<T> *dataField;				/**< A scalar field containing field data (PLAN TO REMOVE) */
	ITL_field_regular<int>* binData;		/**< A scalar field containing histogram bins corresponding to field points. (created or deleted elsewhere) */

	T histogramMin;							/**< Lower limit of histogram range (PLAN TO REMOVE) */
	T histogramMax;							/**< Upper limit of histogram range (PLAN TO REMOVE) */
	bool histogramRangeSet;					/**< Flag indicating the status of histogram range (PLAN TO REMOVE) */

	float globalEntropy;					/**< Value of computed global entropy of the field. */
	int *freqList;							/**< Frequency list (non-normalized) for the entire field used in entropy computation. */
	float* normFreqList;   					/**< Frequency list (normalized) for the entire field used in entropy computation. */
	
	// ADD-BY-ABON 11/07/2011
	ITL_histogram *histogram;				/**< Pointer to the histogram class (created or deleted elsewhere). */
	int nBin;								/**< Number of bins used for histogram computation. */

public:

	/**
	 * Old style constructor.(PLAN TO REMOVE)
	 */
	ITL_globalentropy( ITL_field<T> *f, ITL_histogram *hist )
	{
		dataField = f;
		binData = NULL;
		histogramRangeSet = false;
		histogram = hist;

		freqList = NULL;
		normFreqList = NULL;

	}// End constructor

	/**
	 * Constructor
	 */
	ITL_globalentropy( ITL_field_regular<int> *binF, ITL_histogram *hist, int nbin )
	{
		binData = binF;
		histogram = hist;
		nBin = nbin;

		freqList = NULL;
		normFreqList = NULL;

		freqList = new int[nBin];
		normFreqList = new float[nBin];

	}// End constructor

	/**
	 * Destructor
	 */
	~ITL_globalentropy()
	{
		if( freqList != NULL ) delete [] freqList;
		if( normFreqList != NULL ) delete [] normFreqList;
	}

	/**
	 * Global Entropy computation function.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBin Number of bins used in histogram computation/Number of sample points in KDE estimation.
	 * @paran toNormalize boolean flag, TRUE indicates that the computed entropy value needs to be normalized
	 */
	void
	computeGlobalEntropyOfField( bool toNormalize, int method = 0 )
	{
		if( method == 0 )
		{
			// Check if frequencies are already computed from bin Ids
			// If not, only then compute frequencies from bin Ids
			computeHistogramFrequencies();
			printf( "hi2\n" );
			
			// Compute entropy from frequency list			
			globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( freqList,
										  binData->getSize(), nBin, toNormalize );
	
		}		
		else
		{
			//globalEntropy = ITL_entropycore::computeEntropy_KDEBased(
			//								 binData->getDataFull(),
			//								 binData->getSize(),
			//								 0.0, toNormalize );
		}
	}// End function

	/**
	 * Global entropy computation function for unstructured grid.
	 * Creates a scalar field that contains entropy at each grid vertex.
	 * @param nBins Number of bins used in histogram computation/Number of sample points in KDE estimation.
	 * @paran toNormalize boolean flag, TRUE indicates that the computed entropy value needs to be normalized
	 */
	void
	computeGlobalEntropyOfField_Unstructured( int nBin, bool toNormalize)
	{
		ITL_grid_tetrahedral<SCALAR>* uGrid = dynamic_cast<ITL_grid_tetrahedral<SCALAR>*>(dataField->getGrid());

		int i=0,j=0;

		double rangemin;
		double rangemax;
		rangemin = rangemax = uGrid->vertexList[0].f;
		for(i = 1; i < uGrid->getSize(); i++)
		{
			rangemax = rangemax > uGrid->vertexList[i].f ? rangemax : uGrid->vertexList[i].f;
			rangemin = rangemin < uGrid->vertexList[i].f ? rangemin : uGrid->vertexList[i].f;
		}

		double binWidth = (rangemax - rangemin) / (float)nBin;
		histogramMin = rangemin;
		histogramMax = rangemax;
			
		float* hist = new float[nBin];
		memset(hist, 0, sizeof(float) * nBin);

		//process cell by cell: calculate the contribution of cell i to bin j
		double func[21]= {0};
		double min=0,max=0;
		int posmin=0,posmax=0;
		double tt = 0;
		ITL_grid_tetrahedral<SCALAR>* tetGrid = dynamic_cast<ITL_grid_tetrahedral<SCALAR>*>(uGrid);
		for (int i = 0; i < uGrid->nCell; ++i)
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

		globalEntropy = ITL_entropycore::computeEntropy_HistogramBased( hist, nBin, toNormalize );
	}// End function

	/**
	 * Histogram computation from the bin IDs.
	 */
	void
	computeHistogramFrequencies()
	{
		assert( binData != NULL );
		assert( freqList != NULL );
		assert( normFreqList != NULL );

		printf( "hi 1\n" );
		ITL_histogrammapper<int>::computeHistogramFrequencies( &binData, freqList, nBin );
		printf( "hi 2\n" );

		// Normalize the frequencies
		float nPoint = (float)binData->getSize();
		for( int i=0; i<nBin; i++ )
			normFreqList[i] = freqList[i] / nPoint;

	}

	/**
	 * Histogram frequency data accessor function.
	 * Returns pointer to integer array storing the histogram freqencies.
	 */
	void
	getHistogramFrequencies( int *freqlist )
	{
		assert( freqList != NULL );
		assert( freqlist != NULL );

		// Copy frequency for each bin
		memcpy( freqlist, this->freqList, sizeof( int ) * nBin );	

	}// end function

	/**
	 * Histogram frequency data accessor function.
	 * Returns pointer to integer array storing the histogram freqencies.
	 */
	void getHistogramFrequencies( float *normfreqlist )
	{
		assert( freqList != NULL );
		assert( normfreqlist != NULL );

		// Copy frequency for each bin
		memcpy( normfreqlist, normFreqList, sizeof( float ) * nBin );

	}// end function

	/**
	 * Entropy field accessor function.
	 * Returns computed entropy field.
	 * @return pointer to entropy field.
	 */
	float getGlobalEntropy()
	{
		return this->globalEntropy;
	}// end function
};

#endif
/* ITL_GLOBALENTROPY_H_ */
