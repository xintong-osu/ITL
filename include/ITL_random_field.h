/* ITL_random_field.h
 *
 *  Created on: Jul 12, 2011
 *      Author: leeten
 */

#ifndef ITL_RANDOM_FIELD_H_
#define ITL_RANDOM_FIELD_H_

#include <vector>
using namespace std;

// ADD-BY-LEETEN 08/06/2011-BEGIN
#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
#include <netcdf.h>
// ADD-BY-LEETEN 08/12/2011-BEGIN
#else	// #ifndef	WITH_PNETCDF
	#include <pnetcdf.h>
#endif	// #ifndef	WITH_PNETCDF
// ADD-BY-LEETEN 08/12/2011-END
// ADD-BY-LEETEN 08/06/2011-END
#include <math.h>

// ADd-BY-LEETEN 07/22/2011-BEGIN
#include "ITL_Range.h"
// ADd-BY-LEETEN 07/22/2011-END

#include "ITL_histogram.h"	// ADD-BY-LEETEN 12/02/2011

// ADD-BY-LEETEN 07/23/2011-BEGIN
#include "ITL_SphereSpace.h"
// ADD-BY-LEETEN 07/23/2011-END

#include "liblog.h"
#include "libbuf.h"
#include "libbuf2d.h"
#include "libbuf3d.h"

// ADD-BY-LEETEN 08/06/2011-BEGIN
#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
	#define ASSERT_NETCDF(nc_stmt)	\
	{	\
		int iNcError;	\
		ASSERT_OR_LOG(	\
				NC_NOERR == (iNcError = nc_stmt),	\
				fprintf(stderr, "NetCDF Error: %s.", nc_strerror(iNcError)));	\
	}

// ADD-BY-LEETEN 08/12/2011-BEGIN
#else	// #ifndef	WITH_PNETCDF
	#define ASSERT_NETCDF(nc_stmt)	\
	{	\
		int iNcError;	\
		ASSERT_OR_LOG(	\
				NC_NOERR == (iNcError = nc_stmt),	\
				fprintf(stderr, "NetCDF Error: %s.", ncmpi_strerror(iNcError)));	\
	}

#endif	// #ifndef	WITH_PNETCDF
// ADD-BY-LEETEN 08/12/2011-END
// ADD-BY-LEETEN 08/06/2011-END

class ITLRandomField
{
	ITL_histogram *pcHistogram;	// ADD-BY-LEETEN 12/02/2011

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
		int iGlobalId;	// ADD-BY-LEETEN 08/12/2011
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

		// ADD-BY-LEETEN 08/06/2011-BEGIN
		void
		_DumpDimGeometry2Nc(
				int iDim,
				int iNcId,
				int iVarId,
#ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
				const size_t puStart[],
				const size_t puCount[]
// ADD-BY-LEETEN 09/01/2011-BEGIN
#else	// #ifndef WITH_PNETCDF
   				const MPI_Offset puStart[],
   				const MPI_Offset puCount[],
   				const MPI_Offset puStride[]
#endif	// ifndef WITH_PNETCDF
// ADD-BY-LEETEN 09/01/2011-END
		)
		{
			ASSERT_OR_LOG(
				CGeometry::TYPE_RECTANGULAR_GRID == cGeometry.iType,
				fprintf(stderr, "%d: Un-supported geometry type.", cGeometry.iType) );

			ASSERT_OR_LOG(
				iDim < MAX_DIM,
				fprintf(stderr, "%d: Invalid dim.", iDim) );

			const double *pdCoord = cGeometry.cContent.pcAxisCoords[iDim].pdData;
			if( NULL == pdCoord )
			{
			}
			else
			{
				TBuffer<double> pdTemp;
				pdTemp.alloc(this->piDimLengths[iDim]);

				int iBase = cGeometry.cContent.pcAxisCoords[iDim].iBase;
				int iStep = cGeometry.cContent.pcAxisCoords[iDim].iStep;
				for(int i = 0; i < piDimLengths[iDim]; i++)
					pdTemp[i] = pdCoord[iBase + i * iStep];

				#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
				// dump the geometry of the given dim.
				ASSERT_NETCDF(nc_put_vara_double(
						iNcId,
						iVarId,
						puStart,
						puCount,
						&pdTemp[0]));

				// ADD-BY-LEETEN 08/12/2011-BEGIN
				#else	// #ifndef	WITH_PNETCDF
				// MOD-BY-LEETEN 09/01/2011-FROM:
				// LOG_ERROR(fprintf(stderr, "PNetCDF is not fully supported yet."))
				// TO:
				ASSERT_NETCDF(ncmpi_put_vars_double_all(
						iNcId,
						iVarId,
						puStart,
						puCount,
						puStride,
						&pdTemp[0]));
				// MOD-BY-LEETEN 09/01/2011-END
				#endif	// #ifndef	WITH_PNETCDF
				// ADD-BY-LEETEN 08/12/2011-END
			}
		}
		// ADD-BY-LEETEN 08/06/2011-END

		CBlock()
		{
			iGlobalId = -1;	// ADD-BY-LEETEN 08/12/2011
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
		// ADD-BY-LEETEN 08/06/2011-BEGIN
		int iVarId;
		char szName[1024];
		// ADD-BY-LEETEN 08/06/2011-END

		CRange cRange;

		// ADD-BY-LEETEN 08/06/2011-BEGIN
		void _SetName(
				const char* szName
				)
		{
			strcpy(this->szName, szName);
		}

		// constructor
		CDataComponent()
		{
			iVarId = 0;
			szName[0] = '\0';
		}
		// ADD-BY-LEETEN 08/06/2011-END
	};
	TBuffer<CDataComponent>	pcDataComponentInfo;
	int iBoundDataComponent;

	struct CRandomVariable
	{
		// ADD-BY-LEETEN 08/06/2011-BEGIN
	  //! Name of this random variable
		char szName[128];

	  //! ID of the created feature vector
	  int iVarId;
		// ADD-BY-LEETEN 08/06/2011-END

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
					dSample = (double)cSphereSpace.IMapVectorToPatch(pdFeatureVector);
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
					dSample = log10(dSample);
			}
			return dSample;
		}

		//! get the default range of the random variable. This method is designed for vector orientation
		void
		_GetDefaultRange
		(
				double& dMin,
				double& dMax
		  ) const
		{
			dMin = -HUGE_VAL;
			dMax = +HUGE_VAL;
			int iFeatureLength = this->piFeatureVector.USize();
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
		}

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

		// ADD-BY-LEETEN 08/06/2011-BEGIN
		void _SetName(const char *szName)
		{
			strcpy(this->szName, szName);
		}

		const char* SZGetName() const
		{
		return szName;
		}
		// ADD-BY-LEETEN 08/06/2011-END

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
			iFeatureMapping = DEFAULT_FEATURE_MAPPING;

			// ADD-BY-LEETEN 07/23/2011-BEGIN
			cSphereSpace._CopyDefaultMapping();
			// ADD-BY-LEETEN 07/23/2011-END

			// ADD-BY-LEETEN 07/31/2011-BEGIN
			iNrOfBins = ITL_histogram::DEFAULT_NR_OF_BINS;
			// ADD-BY-LEETEN 07/31/2011-END

			// ADD-BY-LEETEN 08/06/2011-BEGIN
			szName[0] = '\0';

			iVarId = 0;
			// ADD-BY-LEETEN 08/06/2011-END
		}
	};
	vector<CRandomVariable*> vcRandomVariables;
	int iBoundRandomVariable;
public:
	int iNrOfGlobalBlocks;	// ADD-BY-LEETEN 08/12/2011

	// ADD-BY-LEETEN 08/06/2011-BEGIN
	enum {
		GLOBAL_TIME_STAMP = -1,
	};
	vector<int> viTimeStamps;

	void
	_AddTimeStamp(
			int iTimeStamp
			);

	int
	IGetNrOfTimeStamps(
			)const;

	// ID of the created NetCDF file
	int iNcId;

	//! Rank of the current process
	int iRank;

	enum {
		NC_DIM_X,
		NC_DIM_Y,
		NC_DIM_Z,
		NC_DIM_T,
		NC_DIM_BLOCK,
		NC_DIM_GLOBAL_TIME,
		NR_OF_NC_DIMS
	};

	int piNcDimIds[NR_OF_NC_DIMS];
	int piNcDimVarIds[CBlock::MAX_DIM];	// variables to record the coordinates
	int iNcTimeVarId;					// varialbe for the time stamps
	static const char* pszNcDimNames[NR_OF_NC_DIMS];

	//! path/filename of the NetCDF file
	char szNetCdfPathFilename[1024];

	//! ID of the random variables for block-wise global entropy
	vector<int> vcBlockGlobalEntropyId;
	void
	_AddBlockGlobalEntropy
	(
			const int iRvId
			);

	void
	_CreateNetCdf
	(
			const char *szPath,
			const char *szFilenamePrefix
	);

	void
	_CloseNetCdf
	(
	);

	void
	_DumpBlockGeometry2NetCdf
	(
			int iBlockId
	);

	void
	  _DumpData2NetCdf(
			   );

	void
	  _DumpRandomSamples2NetCdf
	  (
	   const int iRandomVariable
	   );

	const int
	IGetNrOfBlocks
	(
			) const;

	const int
	IGetNrOfDataComponents
	(
			) const;

	const int
	  IGetNrOfRandomVariables(
				  )const;
	// ADD-BY-LEETEN 08/06/2011-END

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

	// ADD-BY-LEETEN 08/12/2011-BEGIN
	void
	_MapBlock2GlobalId
	(
			const int piLocal2GlobalMapping[],
			const bool bIs1Based
			);
	// ADD-BY-LEETEN 08/12/2011-END

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
		const int iFeatureMapping
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
