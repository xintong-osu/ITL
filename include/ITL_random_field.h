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

#include "liblog.h"
#include "libbuf.h"
#include "libbuf2d.h"
#include "libbuf3d.h"

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

	struct CDataComponent
	{
		CRange cRange;
	};
	TBuffer<CDataComponent>	pcDataComponentInfo;
	int iBoundDataComponent;

	struct CRandomVariable
	{
		TBuffer<int> piFeatureVector;
		CRange cRange;
		bool bIsUsingOrientation;

		CRandomVariable()
		{
			bIsUsingOrientation = false;
		}
	};
	vector<CRandomVariable> vcRandomVariables;
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
		const bool bIsUsingOrientation
	);

	void
	_DumpBoundBlockFeatureVector
	(
		const int iRandomVariable,
		const char* szGeomPathFilename
	);

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
