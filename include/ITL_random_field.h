/* ITL_random_field.h
 *
 *  Created on: Jul 12, 2011
 *      Author: leeten
 */

#ifndef ITL_RANDOM_FIELD_H_
#define ITL_RANDOM_FIELD_H_

#include <vector>
using namespace std;

#include <math.h>

// ADd-BY-LEETEN 07/22/2011-BEGIN
#include "ITL_Range.h"
// ADd-BY-LEETEN 07/22/2011-END

// ADD-BY-LEETEN 07/23/2011-BEGIN
#include "ITL_SphereSpace.h"
// ADD-BY-LEETEN 07/23/2011-END

#include "liblog.h"
#include "libbuf.h"
#include "libbuf2d.h"
#include "libbuf3d.h"

#if	0	// DEL-BY-LEETEN 07/22/2011-BEGIN
	struct CRange{
		double dMin;
		double dMax;

		CRange()
		{
			dMin = -HUGE_VAL;
			dMax = +HUGE_VAL;
		}

		void _Set(double dMin, double dMax)
		{
			this->dMin = dMin;
			this->dMax = dMax;
		}

		double DClamp(double d)
		{
			return min(max(d, dMin), dMax);
		}
	};
#endif	// DEL-BY-LEETEN 07/22/2011-END

class ITLRandomField
{
	//! structure for accessing elements passed from the application
	struct CArray
	{
		//! pointer to the beginning of the pool passed by the application
		const double *pdData;

		//! offset (0-based) of the first element in this pool
		int iBase;

		//! step between consecutive elements
		int iStep;

		void
		_Set
		(
			//! pointer to the beginning of the pool passed by the application
			const double *pdData,

			//! offset (0-based) of the first element in this pool
			int iBase,

			//! step between consecutive elements
			int iStep
		)
		{
			this->pdData = pdData;
			this->iBase = iBase;
			this->iStep = iStep;
		}
	};

	TBuffer2D<CArray> p2DcBlockData;

	//! a union to store the geometry of a block
	struct CBlock
	{
		enum {
			MAX_DIM = 4
		};
		int piDimLengths[MAX_DIM];

		struct CGeometry
		{
			enum
			{
				TYPE_RECTANGULAR_GRID,
				TYPE_CURVELINEAR_GRID,
				NR_OF_TYPES
			};

			int iType;
			union CContent{
				CArray pcAxisCoords[MAX_DIM];
			} cContent;

			CGeometry()
			{
				iType = TYPE_RECTANGULAR_GRID;
			}

			void
			_SetDimCoords
			(
				const int iDimId,
				double *pdCoord,
				const int iBase,
				const int iStep
			)
			{
				ASSERT_OR_LOG(
					0 <= iDimId && iDimId < MAX_DIM,
					fprintf(stderr, "%d: invalid dim.", iDimId));
				ASSERT_OR_LOG(
					iType == TYPE_RECTANGULAR_GRID,
					fprintf(stderr, "%d: non-matched block type.", iType));
				cContent.pcAxisCoords[iDimId]._Set(pdCoord, iBase, iStep);
			}
		} cGeometry;

		void _SetBlockSize
		(
			const int iDim,
			const int piDimLengths[]
		)
		{
			for(int d = 0; d < MAX_DIM; d++)
			{
				int iDimLength = ( d < iDim )?piDimLengths[d]:1;
				ASSERT_OR_LOG(
						iDimLength>=1,
						fprintf(stderr, "%d: Invalid length for dim %d", iDimLength, d) );
				this->piDimLengths[d] = iDimLength;
			}
		}

		void
		_DumpGeometry(
			const char* szGeomLogPathFilename
		)
		{
			ASSERT_OR_LOG(
				CGeometry::TYPE_RECTANGULAR_GRID == cGeometry.iType,
				fprintf(stderr, "%d: Un-supported geometry type.", cGeometry.iType) );

			// dump the geometry per axis
			FILE *fpGeomLog = fopen(szGeomLogPathFilename, "wt");
			ASSERT_OR_LOG(NULL != fpGeomLog, perror(szGeomLogPathFilename));

			// dump the coordinate.
			// NOTE: currently only the 1st 3 dimensions are dumped
			for(int d = 0; d < 3; d++)
			{
				fprintf(fpGeomLog, "%d\n", piDimLengths[d]);

				double dMin;
				double dMax;
				double dPrev;

				dMin =+HUGE_VAL;
				dMax =-HUGE_VAL;
				dPrev = 0.0;
				for(int i = 0; i < piDimLengths[d]; i++)
				{
					const double *pdCoord = cGeometry.cContent.pcAxisCoords[d].pdData;
					int iBase = cGeometry.cContent.pcAxisCoords[d].iBase;
					int iStep = cGeometry.cContent.pcAxisCoords[d].iStep;
					double dC = pdCoord[iBase + i * iStep];
					double dD = dC - dPrev;
					dPrev = dC;
					fprintf(fpGeomLog, "%e %e\n", dC, dD);
					dMin = min(dMin, dC);
					dMax = max(dMax, dC);
				}
				fprintf(fpGeomLog, "%e - %e\n", dMin, dMax);
			}
			fclose(fpGeomLog);
		}

		CBlock()
		{
			for(int d = 0; d < MAX_DIM; d++)
				piDimLengths[d] = 1;
		}
	};
	TBuffer<CBlock> pcBlockInfo;
	int iBoundBlock;

	// ADd-BY-LEETEN 07/22/2011-BEGIN
	// make CDataComponent accessible
public:
	// ADd-BY-LEETEN 07/22/2011-END

	struct CDataComponent
	{
		CRange cRange;
	};
	TBuffer<CDataComponent>	pcDataComponentInfo;
	int iBoundDataComponent;

	struct CRandomVariable
	{
		// ADd-BY-LEETEN 07/22/2011-BEGIN
		enum {
			//! The default feature mapping
			DEFAULT_FEATURE_MAPPING,

			//! Use the magnitude of the feature vector
			FEATURE_MAGNITUDE = DEFAULT_FEATURE_MAPPING,

			//! Use the raw value of the feature vector, which is only supported for 1D scalar fields
			FEATURE_RAW_VALUE,

			//! Use the orientation of the feature vector. It is only supported for 2D/3D vector fields
			FEATURE_ORIENTATION,

			//! Use the magnitude scale of the feature vector. Essentially, the feature vector is the logarithm of the magnitude
			FEATURE_MAGNITUDE_SCALE,

			//! # support mappings from feature vector to random variable
			NR_OF_FEATURE_MAPPINGS
		};
		// ADd-BY-LEETEN 07/22/2011-END

	  // ADD-BY-LEETEN 07/23/2011-BEGIN
	  CSphereSpace cSphereSpace;
	  // ADD-BY-LEETEN 07/23/2011-END

		TBuffer<int> piFeatureVector;
		CRange cRange;
		#if	0	// MOD-BY-LEETEN 07/22/2011-FROM:
			//! A flag whether the orientation of the feature vector is used as the random variable, which only can be true for 2D or 3D
			bool bIsUsingOrientation;
		#else	// MOD-BY-LEETEN 07/22/2011-TO:
		int iFeatureMapping;


		static
		int
		IConvertStringToFeatureMapping(const char* szString)
		{
			int iFeatureMapping = DEFAULT_FEATURE_MAPPING;

			if( !strncmp(szString, "abs", 3) )
				iFeatureMapping = FEATURE_MAGNITUDE;
			else
			if( !strncmp(szString, "raw", 3) )
				iFeatureMapping = FEATURE_RAW_VALUE;
			else
			if( !strncmp(szString, "dir", 3) )
				iFeatureMapping = FEATURE_ORIENTATION;
			else
			if( !strncmp(szString, "log", 3) )
				iFeatureMapping = FEATURE_MAGNITUDE_SCALE;
			else
				ASSERT_OR_LOG(false, fprintf(stderr, "%s: Unknown feature mapping.", szString));

			return iFeatureMapping;
		}

		void
		_Set
		(
			const int iFeatureLength,
			const int piFeatureVector[],
			const int iFeatureMapping = DEFAULT_FEATURE_MAPPING
		)
		{
			ASSERT_OR_LOG(
					DEFAULT_FEATURE_MAPPING <= iFeatureMapping && iFeatureMapping < NR_OF_FEATURE_MAPPINGS,
					fprintf(stderr, "%d: unsupported feature mapping.", iFeatureMapping)
			);
			this->iFeatureMapping = iFeatureMapping;
			this->piFeatureVector.alloc(iFeatureLength);
			for(int f = 0; f < iFeatureLength; f++)
				this->piFeatureVector[f] = piFeatureVector[f];
		}

		//! Convert the input feature vector to the random variable
		double
		DGetRandomVariable
		(
			const double pdFeatureVector[]

		) const
		{
			int iFeatureLength = this->piFeatureVector.USize();
			double dSample = 0.0;
			switch(iFeatureLength)
			{
			case 1:
			{
				if( FEATURE_RAW_VALUE == iFeatureMapping )
				{
					dSample = pdFeatureVector[0];
					break;
				}
				// otherwise, enter the following part...
			}

			case 2:
			{
				if( FEATURE_ORIENTATION == iFeatureMapping )
				{
					dSample = atan2(pdFeatureVector[1], pdFeatureVector[0]);
					break;
				}
				// otherwise, enter the following part...
			}

			case 3:
			{
				if( FEATURE_ORIENTATION == iFeatureMapping )
				{
					// TEST only...
				    // MOD-BY-LEETEN 07/23/2011-FROM:
						// dSample = (double)(rand()%360);
				    // TO:
					dSample = (double)cSphereSpace.IMapVectorToPatch(pdFeatureVector);
					// MOD-BY-LEETEN 07/23/2011-END
					break;
				}
				// otherwise, enter the following part...
			}

			default:
				ASSERT_OR_LOG(
					FEATURE_MAGNITUDE == iFeatureMapping ||
					FEATURE_MAGNITUDE_SCALE == iFeatureMapping,
					fprintf(stderr, "%d: unknown feature mapping.", iFeatureMapping));

				dSample = 0.0;
				for(int f = 0; f < iFeatureLength; f++)
				{
					dSample += pdFeatureVector[f] * pdFeatureVector[f];
				}
				dSample = sqrt(dSample);

				if( FEATURE_MAGNITUDE_SCALE == iFeatureMapping )
					// MOD-BY-LEETEN 08/05/2011-FROM:
						// dSample = log(dSample);
					// TO:
					dSample = log10(dSample);
					// MOD-BY-LEETEN 08/05/2011-END
			}
                        #if 0  // DEL-BY-LEETEN 07/23/2011-BEGIN
			if( FEATURE_MAGNITUDE_SCALE == iFeatureMapping )
				dSample = log(dSample);
                        #endif // DEL-BY-LEETEN 07/23/2011-END
			return dSample;
		}

		//! get the default range of the random variable. This method is designed for vector orientation
		void
		_GetDefaultRange
		(
				double& dMin,
				double& dMax
		// MOD-BY-LEETEN 07/31/2011-FROM:
			// )
		// TO:
		  ) const
		// MOD-BY-LEETEN 07/31/2011-END
		{
			dMin = -HUGE_VAL;
			dMax = +HUGE_VAL;
			int iFeatureLength = this->piFeatureVector.USize();
			#if 0 // MOD-BY-LEETEN 07/31/2011-FROM:
				if(FEATURE_ORIENTATION == iFeatureMapping)
				{
					switch(iFeatureLength)
					{
					case 2:	dMin = -M_PI;	dMax = +M_PI;	break;

					// MOD-BY-LEETEN 07/23/2011-FROM:
					// case 3:	dMin = 0.0;		dMax = 360;		break;	// the range is from 0 to the number of patches on the semi-sphere
					// TO:
					case 3: dMin = 0.0; dMax = cSphereSpace.IGetNrOfPatches(); break;
					// MOD-BY-LEETEN 07/23/2011-END
					}
				}
			#else	// MOD-BY-LEETEN 07/31/2011-TO:
			switch(iFeatureMapping)
			  {
			  case FEATURE_ORIENTATION:
				switch(iFeatureLength)
				{
				case 2:	dMin = -M_PI;	dMax = +M_PI;	break;
				case 3: dMin = 0.0; dMax = cSphereSpace.IGetNrOfPatches(); break;
				}
				break;
			  case FEATURE_MAGNITUDE:
			    dMin = 0;
			    dMax = HUGE_VAL;
			    break;
			  }
			#endif	// MOD-BY-LEETEN 07/31/2011-END
		}
		#endif	// MOD-BY-LEETEN 07/22/2011-END

		// ADD-BY-LEETEN 07/31/2011-BEGIN
		int iNrOfBins;

		void
		_SetNrOfBins
		(
			const int iNrOfBins
		)
		{
			this->iNrOfBins = iNrOfBins;
		}

		int
		IGetNrOfBins
		(
		) const
		{
		  return iNrOfBins;
		}

		double
		DMapToBin
		(
				double dSample
		) const
		{
				// obtain the default range of the random variable, which is especially useful for orientation
				double dDefaultMin, dDefaultMax;
				_GetDefaultRange(dDefaultMin, dDefaultMax);
				double dHistMin, dHistMax;
				dHistMin = max(cRange.dMin, dDefaultMin);
				dHistMax = min(cRange.dMax, dDefaultMax);

				double dBin = (dSample - dHistMin)/(dHistMax - dHistMin);
				dBin = min(1.0, max(0.0, dBin));
				dBin *= (double)iNrOfBins;
				return dBin;
		}
		// ADD-BY-LEETEN 07/31/2011-END

		CRandomVariable()
		{
			// MOD-BY-LEETEN 07/22/2011-FROM:
				// bIsUsingOrientation = false;
			// TO:
			iFeatureMapping = DEFAULT_FEATURE_MAPPING;
			// MOD-BY-LEETEN 07/22/2011-END

			// ADD-BY-LEETEN 07/23/2011-BEGIN
			cSphereSpace._CopyDefaultMapping();
			// ADD-BY-LEETEN 07/23/2011-END

			// ADD-BY-LEETEN 07/31/2011-BEGIN
			iNrOfBins = ITL_histogram::DEFAULT_NR_OF_BINS;
			// ADD-BY-LEETEN 07/31/2011-END
		}
	};
	// MOD-BY-LEETEN 07/22/2011-FROM:
		// vector<CRandomVariable> vcRandomVariables;
	// TO:
	vector<CRandomVariable*> vcRandomVariables;
	// MOD-BY-LEETEN 07/22/2011-END
	int iBoundRandomVariable;
public:
	/////////////////////////////////////////////////
	int IGetBoundBlock
	(
	);

	CBlock&
	CGetBlock
	(
		const int iBlock
	);

	const CBlock&
	CGetBlock
	(
		const int iBlock
	) const;


	CBlock&
	CGetBoundBlock
	(
	);

	const CBlock&
	CGetBoundBlock
	(
	) const;

	void
	_BindBlock
	(
		int iBlock
	);

	void
	_DumpGeometry
	(
		const char* szGeomLogPathFilename
	);

	/////////////////////////////////////////////////
	CDataComponent&
	CGetDataComponent
	(
		const int iDataComponent
	);

	const CDataComponent&
	CGetDataComponent
	(
		const int iDataComponent
	) const;

	CDataComponent&
	CGetBoundDataComponent
	(
	);

	const CDataComponent&
	CGetBoundDataComponent
	(
	) const;

	void
	_BindDataComponent
	(
		int iDataComponent
	);

	/////////////////////////////////////////////////
	void
	_SetBoundArray
	(
		const double pdData[],
		const int iBase,
		const int iStep
	);

	CArray&
	CGetArray
	(
		int iBlock,
		int iDataComponent
	);

	const CArray&
	CGetArray
	(
		int iBlock,
		int iDataComponent
	)const;

	CArray&
	CGetBoundArray
	(
	);

	const CArray&
	CGetBoundArray
	(
	) const;

	/////////////////////////////////////////////////
	// management of random variables
	void
	_AddRandomVariable
	(
		int *piRvId
	);

	void
	_CheckBoundRandomVariable
	(
	);

	void
	_BindRandomVariable
	(
		const int iRvId
	);

	CRandomVariable&
	CGetRandomVariable
	(
		const int iRandomVariable
	);

	const CRandomVariable&
	CGetRandomVariable
	(
		const int iRandomVariable
	) const;

	CRandomVariable&
	CGetBoundRandomVariable
	(
	);

	const CRandomVariable&
	CGetBoundRandomVariable
	(
	) const;

	void
	_SetFeatureVector
	(
		const int iFeatureLength,
		const int piFeatureVector[],
		// MOD-BY-LEETEN 07/22/2011-FROM:
			// const bool bIsUsingOrientation
		// TO:
		const int iFeatureMapping
		// MOD-BY-LEETEN 07/22/2011-END
	);

	// ADD-BY-LEETEN 07/22/2011-BEGIN
	//! Search for the domain range across all processes and
	// automatically assign the range for the random variable
	void
	_UseDomainRange
	(
		const int iRandonVariable
	);
	// ADD-BY-LEETEN 07/22/2011-END

	void
	_DumpBoundBlockFeatureVector
	(
		const int iRandomVariable,
		const char* szGeomPathFilename
	);

	// ADD-BY-LEETEN 07/22/2011-BEGIN
	void
	_CollectRandomSamplesInBlock
	(
		const int iBlock,
		const int iRandomVariable,
		float *pfRandomSamples,
		// ADD-BY-LEETEN 07/31/2011-BEGIN
		const bool bIsMappingToBins = false
		// ADD-BY-LEETEN 07/31/2011-END
	);
	// ADD-BY-LEETEN 07/22/2011-END

	///////////////////////////////////////////////////////////////////
	void
	_ComputeEntorpyInBoundBlock
	(
		int iRandomVariable,
		const char *szEntropyLogPathFilename
	);

	void
	_ComputeEntorpyFieldInBoundBlock
	(
		const int iRandomVariable,
		const int iDim,
		const double pdLocalNeighborhood[],
		const char *szEntropyPathPrefix
	);

	// ADD-BY-LEETEN 07/31/2011-BEGIN
    void
	_ComputeJointEntorpyInBoundBlock
	(
			int iRandomVariable1,
			int iRandomVariable2,
			const char *szEntropyLogPathFilename
	);
	// ADD-BY-LEETEN 07/31/2011-END

	///////////////////////////////////////////////////////////////////
	void
	_Create
	(
		int iNrOfBlocks,
		int iNrOfDataComponents
	);

	ITLRandomField();
	virtual ~ITLRandomField();
};

#endif /* ITL_RANDOM_FIELD_H_ */
