/*
 * ITL_random_field.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: leeten
 */

	// ADD-BY-LEETEN 08/06/2011-BEGIN
	// DEL-BY-LEETNE 08/12/2011	#include <netcdf.h>
	// ADD-BY-LEETEN 08/06/2011-END

	#include "ITL_header.h"
	#include "ITL_base.h"
	#include "ITL_ioutil.h"
	#include "ITL_vectormatrix.h"
	#include "ITL_histogramconstants.h"
	#include "ITL_localentropy.h"
	#include "ITL_globalentropy.h"
	#include "ITL_localjointentropy.h"
	#include "ITL_globaljointentropy.h"
	#include "ITL_field_regular.h"

	#include "ITL_random_field.h"

// ADD-BY-LEETEN 08/06/2011-BEGIN

// declaration of the static members
const char*
ITLRandomField::pszNcDimNames[] =
{
	"x", "y", "z", "local_time", "block", "time"
};
// ADD-BY-LEETEN 08/06/2011-END

// ADD-BY-LEETEN 08/12/2011-BEGIN
void
ITLRandomField::_MapBlock2GlobalId
(
		const int piLocal2GlobalMapping[],
		const bool bIs1Based
		)
{
	int iMaxGlobalId = -1;
	for(int b = 0; b < IGetNrOfBlocks(); b++)
	{
		int iGlobalId = piLocal2GlobalMapping[b];
		if( true == bIs1Based )
			iGlobalId--;	// from 1-based to 0-based
		this->CGetBlock(b).iGlobalId = iGlobalId;

		// update the maximal global ID in the local side
		iMaxGlobalId = max(iMaxGlobalId, iGlobalId);
	}

	// from ID to #blocks
	iMaxGlobalId++;

	// collect the max. #global blocks to proc. 0 and then broadcast too all others
#if	0	// MOD-BY-LEETEN 08/29/2011-FROM:
	ASSERT_OR_LOG(MPI_SUCCESS == MPI_Reduce(&iMaxGlobalId, &iNrOfGlobalBlocks, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD), "");
	if( 0 == iRank )
		ASSERT_OR_LOG(MPI_SUCCESS == MPI_Bcast(&iNrOfGlobalBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD), "");
#else   // MOD-BY-LEETEN 08/29/2011-TO:
	ASSERT_OR_LOG(MPI_SUCCESS == MPI_Allreduce(MPI_IN_PLACE, &iNrOfGlobalBlocks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD), "");
#endif  // MOD-BY-LEETEN 08/29/2011-END
}
// ADD-BY-LEETEN 08/12/2011-END

// ADD-BY-LEETEN 07/22/2011-BEGIN
void
ITLRandomField::_UseDomainRange
(
	const int iRandomVariable
)
{
	CRandomVariable& cRandomVariable = CGetRandomVariable(iRandomVariable);
	// DEL-BY-LEETEN 08/06/2011-BEGIN
		// const int iFeatureLength = cRandomVariable.piFeatureVector.USize();
	// DEL-BY-LEETEN 08/06/2011-END

	// scan through all block to find the maximal value
	float fLocalMin = +HUGE_VALF;
	float fLocalMax = -HUGE_VALF;
	for(int b = 0; b < (int)this->pcBlockInfo.USize(); b++)
	{
		const CBlock& cBlock = CGetBlock(b);

		TBuffer3D<float> p3DfSamples;
		p3DfSamples.alloc(
			cBlock.piDimLengths[0],
			cBlock.piDimLengths[1],
			cBlock.piDimLengths[2]
			);

		// compute the #cells
		int iNrOfCells = p3DfSamples.USize();

		_CollectRandomSamplesInBlock(b, iRandomVariable, &p3DfSamples[0]);
		for(int c = 0; c < iNrOfCells; c++)
		{
			float fSample = p3DfSamples[c];
			fLocalMin = min(fLocalMin, fSample);
			fLocalMax = max(fLocalMax, fSample);
		}
	}

	#if	0	// DEL-BY-LEETEN 08/12/2011-BEGIN
		int iRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
	#endif	// DEL-BY-LEETEN 08/12/2011-END

#if	0	// MOD-BY-LEETEN 08/29/2011-FROM:
	float fDomainMax;
	MPI_Reduce(&fLocalMax, &fDomainMax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	float fDomainMin;
	MPI_Reduce(&fLocalMin, &fDomainMin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

	#if 0	// DEL-BY-LEETEN 08/06/2011-BEGIN
		if( 0 == iRank )
		{
			LOG_VAR_TO_ERROR(fDomainMin);
			LOG_VAR_TO_ERROR(fDomainMax);
		}
	#endif	// DEL-BY-LEETEN 08/06/2011-END

	MPI_Bcast(&fDomainMin, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&fDomainMax, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#else   // MOD-BY-LEETEN 08/29/2011-TO:
	float fDomainMin;
	float fDomainMax;
	MPI_Allreduce(MPI_IN_PLACE, &fDomainMin, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &fDomainMax, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#endif  // MOD-BY-LEETEN 08/29/2011-END

	cRandomVariable.cRange._Set((double)fDomainMin, (double)fDomainMax);
}
// ADD-BY-LEETEN 07/22/2011-END

// ADD-BY-LEETEN 08/06/2011-BEGIN
const int 
ITLRandomField::IGetNrOfRandomVariables(
				  )const
{
  return (int)vcRandomVariables.size();
}

int
ITLRandomField::IGetNrOfTimeStamps(
		)const
{
	return (int)viTimeStamps.size();
}

void
ITLRandomField::_AddTimeStamp(
		int iTimeStamp
		)
{
	viTimeStamps.push_back(iTimeStamp);
}

const int
ITLRandomField::IGetNrOfBlocks
(
		) const
{
	return (int)this->pcBlockInfo.USize();
}

const int
ITLRandomField::IGetNrOfDataComponents
(
		) const
{
	return (int)this->pcDataComponentInfo.USize();
}
// ADD-BY-LEETEN 08/06/2011-END

ITLRandomField::CBlock&
ITLRandomField::CGetBlock
(
	const int iBlock
)
{
	return pcBlockInfo[iBlock];
}

const ITLRandomField::CBlock&
ITLRandomField::CGetBlock
(
	const int iBlock
) const
{
	return pcBlockInfo[iBlock];
}

int
ITLRandomField::IGetBoundBlock
(
)
{
	return iBoundBlock;
}

ITLRandomField::CBlock&
ITLRandomField::CGetBoundBlock
(
)
{
	return CGetBlock(iBoundBlock);
}

const ITLRandomField::CBlock&
ITLRandomField::CGetBoundBlock
(
) const
{
	return CGetBlock(iBoundBlock);
}

void
ITLRandomField::_BindBlock
(
	const int iBlock
)
{
	ASSERT_OR_LOG(
		0 <= iBlock && iBlock < (int)this->pcBlockInfo.USize(),
		fprintf(stderr, "%d: invalid block.", iBlock));
	iBoundBlock = iBlock;
}

void
ITLRandomField::_DumpGeometry
(
	const char* szGeomLogPathFilename
)
{
	pcBlockInfo[iBoundBlock]._DumpGeometry(szGeomLogPathFilename);
}

////////////////////////////////////////////////
ITLRandomField::CDataComponent&
ITLRandomField::CGetDataComponent
(
	const int iDataComponent
)
{
	return pcDataComponentInfo[iDataComponent];
}

const
ITLRandomField::CDataComponent&
ITLRandomField::CGetDataComponent
(
	const int iDataComponent
) const
{
	return pcDataComponentInfo[iDataComponent];
}

ITLRandomField::CDataComponent&
ITLRandomField::CGetBoundDataComponent
(
)
{
	return CGetDataComponent(iBoundDataComponent);
}

const
ITLRandomField::CDataComponent&
ITLRandomField::CGetBoundDataComponent
(
) const
{
	return CGetDataComponent(iBoundDataComponent);
}

void
ITLRandomField::
_BindDataComponent
(
	int iDataComponent
)
{
	ASSERT_OR_LOG(
		0 <= iDataComponent && iDataComponent < (int)this->pcDataComponentInfo.USize(),
		fprintf(stderr, "%d: invalid data component.", iDataComponent));
	iBoundDataComponent = iDataComponent;
}

/////////////////////////////////////////////////
void
ITLRandomField::_SetBoundArray
(
	const double pdData[],
	const int iBase,
	const int iStep
)
{
	this->CGetBoundArray()._Set(pdData, iBase, iStep);
}

ITLRandomField::CArray&
ITLRandomField::CGetArray
(
	int iBlock,
	int iDataComponent
)
{
	return this->p2DcBlockData.GetElement(iBlock, iDataComponent);
}

const
ITLRandomField::CArray&
ITLRandomField::CGetArray
(
	int iBlock,
	int iDataComponent
)const
{
	return this->p2DcBlockData.GetElement(iBlock, iDataComponent);
}

ITLRandomField::CArray&
ITLRandomField::CGetBoundArray
(
)
{
	return CGetArray(iBoundBlock, iBoundDataComponent);
}

const
ITLRandomField::CArray&
ITLRandomField::CGetBoundArray
(
)const
{
	return CGetArray(iBoundBlock, iBoundDataComponent);
}

//////////////////////////////////////////////////////////////////////////
void
ITLRandomField::_AddRandomVariable
(
	int *piRvId
)
{
	int iRvId = vcRandomVariables.size();
	// MOD-BY-LEETEN 07/22/2011-FROM:
		// vcRandomVariables.push_back(CRandomVariable());
	// TO:
	vcRandomVariables.push_back(new CRandomVariable());
	// MOD-BY-LEETEN 07/22/2011-END
	*piRvId = iRvId;

}

void
ITLRandomField::_CheckBoundRandomVariable
(
)
{
	ASSERT_OR_LOG(
		0 <= iBoundRandomVariable && iBoundRandomVariable < (int)vcRandomVariables.size(),
		fprintf(stderr, "%d: invalid or non-bound random variable", iBoundRandomVariable));
}

void
ITLRandomField::_BindRandomVariable
(
	const int iRandomVariable
)
{
	iBoundRandomVariable = iRandomVariable;
	_CheckBoundRandomVariable();
}

ITLRandomField::CRandomVariable&
ITLRandomField::CGetRandomVariable
(
	const int iRandomVariable
)
{
	// MOD-BY-LEETEN 07/22/2011-FROM:
		// return vcRandomVariables[iRandomVariable];
	// TO:
	return *vcRandomVariables[iRandomVariable];
	// MOD-BY-LEETEN 07/22/2011-END
}

const ITLRandomField::CRandomVariable&
ITLRandomField::CGetRandomVariable
(
	const int iRandomVariable
) const
{
	// MOD-BY-LEETEN 07/22/2011-FROM:
		// return vcRandomVariables[iRandomVariable];
	// TO:
	return *vcRandomVariables[iRandomVariable];
	// MOD-BY-LEETEN 07/22/2011-END
}

ITLRandomField::CRandomVariable&
ITLRandomField::CGetBoundRandomVariable
(
)
{
	return CGetRandomVariable(iBoundRandomVariable);
}

const ITLRandomField::CRandomVariable&
ITLRandomField::CGetBoundRandomVariable
(
) const
{
	return CGetRandomVariable(iBoundRandomVariable);
}


void
ITLRandomField::_SetFeatureVector
(
	const int iFeatureLength,
	const int piFeatureVector[],
	// MOD-BY-LEETEN 07/22/2011-FROM:
		// const bool bIsUsingOrientation
	// TO:
	const int iFeatureMapping
	// MOD-BY-LEETEN 07/22/2011-END
)
{
#if	0	// MOD-BY-LEETEN 07/21/2011-FROM:
	CRandomVariable& cBoundRandomVariable = CGetBoundRandomVariable();
	cBoundRandomVariable.bIsUsingOrientation = bIsUsingOrientation;
	cBoundRandomVariable.piFeatureVector.alloc(iFeatureLength);
	for(int f = 0; f < iFeatureLength; f++)
		cBoundRandomVariable.piFeatureVector[f] = piFeatureVector[f];
#else
	CGetBoundRandomVariable()._Set(iFeatureLength, piFeatureVector, iFeatureMapping);
#endif
}

void
ITLRandomField::_DumpBoundBlockFeatureVector
(
	const int iRandomVariable,
	const char* szGeomPathFilename
)
{
	const int iBoundBlock = IGetBoundBlock();
	const CBlock& cBoundBlock = CGetBoundBlock();

	CRandomVariable& cRandomVariable = CGetRandomVariable(iRandomVariable);
	const int iFeatureLength = cRandomVariable.piFeatureVector.USize();

	// compute the #cells
	int iNrOfCells = 1;
	for(int d = 0;	d < CBlock::MAX_DIM; d++)
	{
		iNrOfCells *= cBoundBlock.piDimLengths[d];
	}

	switch(iFeatureLength)
	{
	case 1:
	{
		TBuffer3D<float> p3DfFeatureScalars;

		p3DfFeatureScalars.alloc(
			cBoundBlock.piDimLengths[0],
			cBoundBlock.piDimLengths[1],
			cBoundBlock.piDimLengths[2]
			);

		for(int f = 0; f < iFeatureLength; f++)
		{
			int iDataComponent = cRandomVariable.piFeatureVector[f];
			CDataComponent& cDataComponent = this->CGetDataComponent(iDataComponent);

			const CArray &cArray = this->CGetArray(iBoundBlock, iDataComponent);
			const double *pdData = cArray.pdData;
			int iBase = cArray.iBase;
			int iStep = cArray.iStep;
			for(int c = 0;	c < iNrOfCells; c++)
			{
				double dValue = (float)pdData[iBase + c * iStep];
				p3DfFeatureScalars[c] = cDataComponent.cRange.DClamp(dValue);
			}
		}
		p3DfFeatureScalars._Save(szGeomPathFilename);
	} break;

	// ADD-BY-LEETEN 07/22/2011-BEGIN
	case 2:
	// ADD-BY-LEETEN 07/22/2011-END
	case 3:
	{
		TBuffer3D<VECTOR3> p3DfFeatureVectors;

		p3DfFeatureVectors.alloc(
			cBoundBlock.piDimLengths[0],
			cBoundBlock.piDimLengths[1],
			cBoundBlock.piDimLengths[2]
			);

		for(int f = 0; f < iFeatureLength; f++)
		{
			int iDataComponent = cRandomVariable.piFeatureVector[f];
			CDataComponent& cDataComponent = this->CGetDataComponent(iDataComponent);

			const CArray &cArray = this->CGetArray(iBoundBlock, iDataComponent);
			const double *pdData = cArray.pdData;
			int iBase = cArray.iBase;
			int iStep = cArray.iStep;
			for(int c = 0;	c < iNrOfCells; c++)
			{
				double dValue = pdData[iBase + c * iStep];
				p3DfFeatureVectors[c][f] = (float)cDataComponent.cRange.DClamp(dValue);
			}
		}
		// ADD-BY-LEETEN 07/22/2011-BEGIN
		if( 2 == iFeatureLength )
			for(int c = 0;	c < iNrOfCells; c++)
				p3DfFeatureVectors[c][2] = 0.0f;
		// ADD-BY-LEETEN 07/22/2011-END

		p3DfFeatureVectors._Save(szGeomPathFilename);
	} break;
	}
}

// ADD-BY-LEETEN 07/22/2011-BEGIN
void
ITLRandomField::_CollectRandomSamplesInBlock
(
	const int iBlock,
	const int iRandomVariable,	// CRandomVariable& cRandomVariable,
	float pfRandomSamples[],
	// ADD-BY-LEETEN 07/31/2011-BEGIN
	const bool bIsMappingToBins
	// ADD-BY-LEETEN 07/31/2011-END
)
{
	const CRandomVariable& cRandomVariable = this->CGetRandomVariable(iRandomVariable);
	int iFeatureLength = cRandomVariable.piFeatureVector.USize();

	// allocate the feature vector
	TBuffer<double> pdFeatureVector;
	pdFeatureVector.alloc(iFeatureLength);

	const CBlock& cBlock = this->CGetBlock(iBlock);
	int iNrOfCells = 1;
	for(int d = 0; d < CBlock::MAX_DIM; d++)
		iNrOfCells *= cBlock.piDimLengths[d];

	for(int c = 0; c < iNrOfCells; c++)
	{
		// collect the feature vector
		for(int f = 0; f < iFeatureLength; f++)
		{
			int iDataComponent = cRandomVariable.piFeatureVector[f];
			const CDataComponent& cDataComponent = this->CGetDataComponent(iDataComponent);

			const CArray &cArray = this->CGetArray(iBlock, iDataComponent);
			const double *pdData = cArray.pdData;
			int iBase = cArray.iBase;
			int iStep = cArray.iStep;
			double dValue = pdData[iBase + c * iStep];
			pdFeatureVector[f] = cDataComponent.cRange.DClamp(dValue);
		}
		double dSample = cRandomVariable.DGetRandomVariable(&pdFeatureVector[0]);
		// MOD-BY-LEETEN 07/31/2011-FROM:
			// pfRandomSamples[c] = (float)cRandomVariable.cRange.DClamp(dSample);
		// TO:
		dSample = cRandomVariable.cRange.DClamp(dSample);
		if( bIsMappingToBins )
		  {
		    dSample = cRandomVariable.DMapToBin(dSample);
		  }
		pfRandomSamples[c] = (float)dSample;
		// MOD-BY-LEETEN 07/31/2011-END
	}
}
// ADD-BY-LEETEN 07/22/2011-END


/////////////////////////////////////////////////////////////////
void
ITLRandomField::_Create
(
	int iNrOfBlocks,
	int iNrOfDataComponents
)
{
	ASSERT_OR_LOG(iNrOfBlocks, fprintf(stderr, "%d: invalide #blocks.", iNrOfBlocks));
	this->pcBlockInfo.alloc(iNrOfBlocks);

	ASSERT_OR_LOG(iNrOfDataComponents, fprintf(stderr, "%d: invalide #data-components.", iNrOfDataComponents));
	// MOD-BY-LEETEN 07/18/2011-FROM:
		// this->pcDataComponentInfo.alloc(iNrOfDataComponents);
	// TO:
	this->pcDataComponentInfo.New(iNrOfDataComponents);
	// MOD-BY-LEETEN 07/18/2011-END

	this->p2DcBlockData.alloc(iNrOfBlocks, iNrOfDataComponents);

	// ADD-BY-LEETEN 08/06/2011-BEGIN
	// get the current rank
	MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
	// ADD-BY-LEETEN 08/06/2011-END
}

// ADD-BY-LEETEN 08/06/2011-BEGIN
/////////////////////////////////////////////////////////////////
void
ITLRandomField::_CreateNetCdf
(
		const char *szPath,
		const char *szFilenamePrefix
)
{
	char szNetCdfPathFilename[1024];
	#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
	sprintf(szNetCdfPathFilename, "%s/%s.rank_%d.nc", szPath, szFilenamePrefix, iRank);

    // Create the file.
    ASSERT_NETCDF(nc_create(
    		szNetCdfPathFilename,
    		NC_CLOBBER,
    		&iNcId));

    // ADD-BY-LEETEN 08/12/2011-BEGIN
	#else	// #ifndef	WITH_PNETCDF
	sprintf(szNetCdfPathFilename, "%s/%s.nc", szPath, szFilenamePrefix);
    ASSERT_NETCDF(ncmpi_create(
    		MPI_COMM_WORLD,
    		szNetCdfPathFilename,
    		NC_CLOBBER,
    		MPI_INFO_NULL,
    		&iNcId));
	#endif	// #ifndef	WITH_PNETCDF
    // ADD-BY-LEETEN 08/12/2011-END

    // find the maximal block dim
    int piBlockDimMaxLengths[CBlock::MAX_DIM];
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		piBlockDimMaxLengths[d] = 0;
	}
    for(int b = 0; b < (int)pcBlockInfo.USize(); b++)
    	for(int d = 0; d < CBlock::MAX_DIM; d++)
    		piBlockDimMaxLengths[d] = max(piBlockDimMaxLengths[d], pcBlockInfo[b].piDimLengths[d]);

    // ADD-BY-LEETEN 08/12/2011-BEGIN
	// collect the max. length of all dim.
#if	0	// MOD-BY-LEETEN 08/29/2011-FROM:
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		int iTemp = piBlockDimMaxLengths[d];
		ASSERT_OR_LOG(MPI_SUCCESS == MPI_Reduce(&iTemp, &piBlockDimMaxLengths[d], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD), "");
		if( 0 == iRank )
			ASSERT_OR_LOG(MPI_SUCCESS == MPI_Bcast(&piBlockDimMaxLengths[d], 1, MPI_INT, 0, MPI_COMM_WORLD), "");
	}
#else   // MOD-BY-LEETEN 08/29/2011-TO:
	ASSERT_OR_LOG(MPI_SUCCESS == MPI_Allreduce(MPI_IN_PLACE, &piBlockDimMaxLengths[0], CBlock::MAX_DIM, MPI_INT, MPI_MAX, MPI_COMM_WORLD), "");
#endif  // MOD-BY-LEETEN 08/29/2011-END

	#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
	for(int d = 0; d < CBlock::MAX_DIM; d++)
		ASSERT_NETCDF(nc_def_dim(
						iNcId,
						pszNcDimNames[d],
						piBlockDimMaxLengths[d],
						&piNcDimIds[d]));

    // Define the block dimension
	ASSERT_NETCDF(nc_def_dim(
    				iNcId,
    				pszNcDimNames[NC_DIM_BLOCK],
    				IGetNrOfBlocks(),
    				&piNcDimIds[NC_DIM_BLOCK]));

	// Define the time dimension with unlimited length
    ASSERT_NETCDF(nc_def_dim(
    				iNcId,
    				pszNcDimNames[NC_DIM_GLOBAL_TIME],
    				NC_UNLIMITED,
    				&piNcDimIds[NC_DIM_GLOBAL_TIME]));

    // Define the variables for the coordinates
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		int piBlockDims[3];
		piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
		piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
		piBlockDims[2] = piNcDimIds[d];

		ASSERT_NETCDF(nc_def_var(
						iNcId,
						pszNcDimNames[d],
						NC_DOUBLE,
						sizeof(piBlockDims) / sizeof(piBlockDims[0]),
						piBlockDims,
						&piNcDimVarIds[d]));
	}

	// define the variable for time stamp
	ASSERT_NETCDF(nc_def_var(
			iNcId,
			pszNcDimNames[NC_DIM_GLOBAL_TIME],
			NC_INT,
			1,
			&piNcDimIds[NC_DIM_GLOBAL_TIME],
	     	&iNcTimeVarId));

	// define the variable for all data components
	for(int c = 0; c < this->IGetNrOfDataComponents(); c++)
	{
		int piBlockDims[6];
		piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
		piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
		piBlockDims[2] = piNcDimIds[NC_DIM_T];
		piBlockDims[3] = piNcDimIds[NC_DIM_Z];
		piBlockDims[4] = piNcDimIds[NC_DIM_Y];
		piBlockDims[5] = piNcDimIds[NC_DIM_X];

		ASSERT_NETCDF(nc_def_var(
						iNcId,
						this->CGetDataComponent(c).szName,
						NC_DOUBLE,
						sizeof(piBlockDims) / sizeof(piBlockDims[0]),
						piBlockDims,
						&this->CGetDataComponent(c).iVarId));
	}

	// define the varaibles for all random variables
	for(int r = 0; r < this->IGetNrOfRandomVariables(); r++)
	{
	  CRandomVariable& cRv = this->CGetRandomVariable(r);
	  int piBlockDims[6];
	  piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
	  piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
	  piBlockDims[2] = piNcDimIds[NC_DIM_T];
	  piBlockDims[3] = piNcDimIds[NC_DIM_Z];
	  piBlockDims[4] = piNcDimIds[NC_DIM_Y];
	  piBlockDims[5] = piNcDimIds[NC_DIM_X];

	  ASSERT_NETCDF(nc_def_var(
				   iNcId,
				   cRv.szName,
				   NC_FLOAT,
				   sizeof(piBlockDims) / sizeof(piBlockDims[0]),
				   piBlockDims,
				   &cRv.iVarId));
	}

	// finish the definition mode
	ASSERT_NETCDF(nc_enddef(
			iNcId));

	// ADD-BY-LEETEN 08/12/2011-BEGIN
	#else	// #ifndef	WITH_PNETCDF
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		ASSERT_NETCDF(ncmpi_def_dim(
						iNcId,
						pszNcDimNames[d],
						piBlockDimMaxLengths[d],
						&piNcDimIds[d]));
	}

    // Define the block dimension
	ASSERT_NETCDF(ncmpi_def_dim(
    				iNcId,
    				pszNcDimNames[NC_DIM_BLOCK],
				// MOD-BY-LEETEN 08/30/2011-FROM:
    				  // IGetNrOfBlocks(),
				// TO:
				iNrOfGlobalBlocks,
				// MOD-BY-LEETEN 08/30/2011-END
    				&piNcDimIds[NC_DIM_BLOCK]));

    // Define the time dimension with unlimited length
    ASSERT_NETCDF(ncmpi_def_dim(
    				iNcId,
    				pszNcDimNames[NC_DIM_GLOBAL_TIME],
    				NC_UNLIMITED,
    				&piNcDimIds[NC_DIM_GLOBAL_TIME]));

    // Define the variables for the coordinates
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		int piBlockDims[3];
		piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
		piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
		piBlockDims[2] = piNcDimIds[d];

		ASSERT_NETCDF(ncmpi_def_var(
						iNcId,
						pszNcDimNames[d],
						NC_DOUBLE,
						sizeof(piBlockDims) / sizeof(piBlockDims[0]),
						piBlockDims,
						&piNcDimVarIds[d]));
	}

	// define the variable for time stamp
	ASSERT_NETCDF(ncmpi_def_var(
			iNcId,
			pszNcDimNames[NC_DIM_GLOBAL_TIME],
			NC_INT,
			1,
			&piNcDimIds[NC_DIM_GLOBAL_TIME],
	     	&iNcTimeVarId));

	// define the variable for all data components
	for(int c = 0; c < this->IGetNrOfDataComponents(); c++)
	{
		int piBlockDims[6];
		piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
		piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
		piBlockDims[2] = piNcDimIds[NC_DIM_T];
		piBlockDims[3] = piNcDimIds[NC_DIM_Z];
		piBlockDims[4] = piNcDimIds[NC_DIM_Y];
		piBlockDims[5] = piNcDimIds[NC_DIM_X];

		ASSERT_NETCDF(ncmpi_def_var(
						iNcId,
						this->CGetDataComponent(c).szName,
						NC_DOUBLE,
						sizeof(piBlockDims) / sizeof(piBlockDims[0]),
						piBlockDims,
						&this->CGetDataComponent(c).iVarId));
	}

	// define the varaibles for all random variables
	for(int r = 0; r < this->IGetNrOfRandomVariables(); r++)
	{
	  CRandomVariable& cRv = this->CGetRandomVariable(r);
	  int piBlockDims[6];
	  piBlockDims[0] = piNcDimIds[NC_DIM_GLOBAL_TIME];
	  piBlockDims[1] = piNcDimIds[NC_DIM_BLOCK];
	  piBlockDims[2] = piNcDimIds[NC_DIM_T];
	  piBlockDims[3] = piNcDimIds[NC_DIM_Z];
	  piBlockDims[4] = piNcDimIds[NC_DIM_Y];
	  piBlockDims[5] = piNcDimIds[NC_DIM_X];

	  ASSERT_NETCDF(ncmpi_def_var(
				   iNcId,
				   cRv.szName,
				   NC_FLOAT,
				   sizeof(piBlockDims) / sizeof(piBlockDims[0]),
				   piBlockDims,
				   &cRv.iVarId));
	}

	// finish the definition mode
	ASSERT_NETCDF(ncmpi_enddef(
			iNcId));
	#endif	// #ifndef	WITH_PNETCDF
	// ADD-BY-LEETEN 08/12/2011-END
	// enter the data mode...
}

void
ITLRandomField::_CloseNetCdf
(
)
{
	if( iNcId > 0 )
	{
		// write the time stamp
		TBuffer<int> piTemp;
		piTemp.alloc(this->IGetNrOfTimeStamps());
		for(int t = 0; t < (int)piTemp.USize(); t++)
			piTemp[t] = this->viTimeStamps[t];

		#ifndef	WITH_PNETCDF		// ADD-BY-LEETEN 08/12/2011
		size_t uStart = 0;
		size_t uCount = piTemp.USize();
		ASSERT_NETCDF(nc_put_vara_int(
				iNcId,
				iNcTimeVarId,
				&uStart,
				&uCount,
				&piTemp[0]));
        /* Close the file. */
	    ASSERT_NETCDF(nc_close(iNcId));

		// ADD-BY-LEETEN 08/12/2011-BEGIN
		#else	// #ifndef	WITH_PNETCDF
	    MPI_Offset uStart = 0;
	    MPI_Offset uCount = piTemp.USize();

		ASSERT_NETCDF(ncmpi_begin_indep_data(iNcId));
	    if( 0 == iRank )
			ASSERT_NETCDF(ncmpi_put_vara_int(
					iNcId,
					iNcTimeVarId,
					&uStart,
					&uCount,
					&piTemp[0]));
		ASSERT_NETCDF(ncmpi_end_indep_data(iNcId));

        /* Close the file. */
	    ASSERT_NETCDF(ncmpi_close(iNcId));
		#endif	// #ifndef	WITH_PNETCDF
		// ADD-BY-LEETEN 08/12/2011-END

	    iNcId = 0;
	}
};

void
ITLRandomField::_DumpBlockGeometry2NetCdf
(
		int iBlockId
		)
{
    /* Write the coordinate variable data. This will put the latitudes
       and longitudes of our data grid into the netCDF file. */
	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		size_t puStart[3];
		size_t puCount[3];

		// time
		puStart[0] = this->IGetNrOfTimeStamps() - 1;
		puCount[0] = 1;

		// block
		puStart[1] = (size_t)iBlockId;
		puCount[1] = 1;

		// spatial
		puStart[2] = 0;
		puCount[2] = (size_t)this->CGetBlock(iBlockId).piDimLengths[d];

		this->CGetBlock(iBlockId)._DumpDimGeometry2Nc(
						d,
						iNcId,
						piNcDimVarIds[d],
						puStart,
						puCount);
	}
}

void
ITLRandomField::_DumpData2NetCdf
(
		)
{
  for(int b = 0; b < IGetNrOfBlocks(); b++)
    {
      const CBlock& cBlock = this->CGetBlock(b);
      int iNrOfCells = 1;
      for(int d = 0; d < CBlock::MAX_DIM; d++)
	iNrOfCells *= cBlock.piDimLengths[d];

      TBuffer<double> pdTemp;
      pdTemp.alloc(iNrOfCells);

      size_t puStart[6];
      size_t puCount[6];
      // time
      puStart[0] = this->IGetNrOfTimeStamps() - 1;
      puCount[0] = 1;

      // block
      puStart[1] = (size_t)b;
      puCount[1] = 1;

      for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
	  puStart[2+d] = 0;
	  puCount[2+d] = (size_t)cBlock.piDimLengths[CBlock::MAX_DIM - 1 - d];
	}

      for(int c = 0; c < IGetNrOfDataComponents(); c++)
	{ 
	  const CDataComponent& cDataComponent = CGetDataComponent(c);
	  const CArray &cArray = this->CGetArray(b, c);
	  const double *pdData = cArray.pdData;
	  int iBase = cArray.iBase;
	  int iStep = cArray.iStep;

	  for(int v = 0; v < iNrOfCells; v++)
	    {
	      double dValue = pdData[iBase + v * iStep];
	      pdTemp[v] = cDataComponent.cRange.DClamp(dValue);
	    }

		#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
	  // dump the geometry of the given dim.
	  ASSERT_NETCDF(nc_put_vara_double(
					   iNcId,
					   cDataComponent.iVarId,
					   puStart,
					   puCount,
					   &pdTemp[0]));

		// ADD-BY-LEETEN 08/12/2011-BEGIN
		#else	// #ifndef	WITH_PNETCDF
		LOG_ERROR(fprintf(stderr, "PNetCDF is not fully supported yet."))
		#endif	// #ifndef	WITH_PNETCDF
		// ADD-BY-LEETEN 08/12/2011-END
	}
    }
}

void
ITLRandomField::_DumpRandomSamples2NetCdf
(
 const int iRandomVariable
		)
{
  for(int b = 0; b < IGetNrOfBlocks(); b++)
    {
      const CBlock& cBlock = this->CGetBlock(b);
      int iNrOfCells = 1;
      for(int d = 0; d < CBlock::MAX_DIM; d++)
	iNrOfCells *= cBlock.piDimLengths[d];

      TBuffer<float> pfTemp;
      pfTemp.alloc(iNrOfCells);

      size_t puStart[6];
      size_t puCount[6];

      // time
      puStart[0] = this->IGetNrOfTimeStamps() - 1;
      puCount[0] = 1;

      // block
      puStart[1] = (size_t)b;
      puCount[1] = 1;

      for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
	  puStart[2+d] = 0;
	  puCount[2+d] = (size_t)cBlock.piDimLengths[CBlock::MAX_DIM - 1 - d];
	}

      _CollectRandomSamplesInBlock(b, iRandomVariable, &pfTemp[0], true);


		#ifndef	WITH_PNETCDF	// ADD-BY-LEETEN 08/12/2011
      ASSERT_NETCDF(nc_put_vara_float(
				      iNcId,
				      CGetRandomVariable(iRandomVariable).iVarId,
				      puStart,
				      puCount,
				      &pfTemp[0]));
		// ADD-BY-LEETEN 08/12/2011-BEGIN
		#else	// #ifndef	WITH_PNETCDF
		LOG_ERROR(fprintf(stderr, "PNetCDF is not fully supported yet."))
		#endif	// #ifndef	WITH_PNETCDF
		// ADD-BY-LEETEN 08/12/2011-END
    }
}

// ADD-BY-LEETEN 08/06/2011-END

//////////////////////////////////////////////////////////////////////////
void
ITLRandomField::_ComputeEntorpyInBoundBlock
(
	int iRandomVariable,
	const char *szEntropyLogPathFilename
)
{
	const int iBoundBlock = IGetBoundBlock();
	const CBlock& cBoundBlock = CGetBoundBlock();

	CRandomVariable& cRandomVariable = CGetRandomVariable(iRandomVariable);
	// DEL-BY-LEETEN 08/06/2011-BEGIN
		// const int iFeatureLength = cRandomVariable.piFeatureVector.USize();
	// DEL-BY-LEETEN 08/06/2011-END

	// compute the #cells
	int iNrOfCells = 1;
	for(int d = 0;	d < CBlock::MAX_DIM; d++)
	{
		iNrOfCells *= cBoundBlock.piDimLengths[d];
	}

	// other setting of the blocks.
	float pfBlockDimLow[CBlock::MAX_DIM];
	float pfBlockDimUp[CBlock::MAX_DIM];
	int	piLowPad[CBlock::MAX_DIM];
	int piHighPad[CBlock::MAX_DIM];
	int piNeighborhood[CBlock::MAX_DIM];
	// DEL-BY-LEETEN 07/31/2011-BEGIN
		// int iDefaultNrOfBins = 360;
	// DEL-BY-LEETEN 07/31/2011-END

	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		pfBlockDimLow[d] = 0.0f;
		pfBlockDimUp[d] = (float)cBoundBlock.piDimLengths[d] - 1;
		piLowPad[d] = 0;
		piHighPad[d] = 0;
		piNeighborhood[d] = 0;
	}

	#if	0	// MOD-BY-LEETEN 07/22/2011-FROM:
		float fEntropy;
		switch(iFeatureLength)
		{
		case 1:
		{
			TBuffer3D<float> p3DfFeatureScalars;

			p3DfFeatureScalars.alloc(
				cBoundBlock.piDimLengths[0],
				cBoundBlock.piDimLengths[1],
				cBoundBlock.piDimLengths[2]
				);

			for(int f = 0; f < iFeatureLength; f++)
			{
				int iDataComponent = cRandomVariable.piFeatureVector[f];
				CDataComponent& cDataComponent = this->CGetDataComponent(iDataComponent);

				const CArray &cArray = this->CGetArray(iBoundBlock, iDataComponent);
				const double *pdData = cArray.pdData;
				int iBase = cArray.iBase;
				int iStep = cArray.iStep;
				for(int c = 0;	c < iNrOfCells; c++)
				{
					double dValue = (float)pdData[iBase + c * iStep];
					p3DfFeatureScalars[c] = cDataComponent.cRange.DClamp(dValue);
				}
			}

			ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
				&p3DfFeatureScalars[0],
				3,
				pfBlockDimLow,
				pfBlockDimUp,
				piLowPad,
				piHighPad,
				piNeighborhood );

			ITL_globalentropy<SCALAR> *globalEntropyComputerForScalar = new ITL_globalentropy<SCALAR>( scalarField );

			// convert each vector into bin index
			globalEntropyComputerForScalar->computeHistogramBinField("scalar", iDefaultNrOfBins);

			// compute the entropy
			globalEntropyComputerForScalar->computeGlobalEntropyOfField(iDefaultNrOfBins, false);

			// save the block-wise entropy
			fEntropy = globalEntropyComputerForScalar->getGlobalEntropy();

			delete globalEntropyComputerForScalar;
		} break;

		case 3:
		{
			TBuffer3D<VECTOR3> p3DfFeatureVectors;

			p3DfFeatureVectors.alloc(
				cBoundBlock.piDimLengths[0],
				cBoundBlock.piDimLengths[1],
				cBoundBlock.piDimLengths[2]
				);

			for(int f = 0; f < iFeatureLength; f++)
			{
				const CArray &cArray = this->CGetArray(iBoundBlock, cRandomVariable.piFeatureVector[f]);
				const double *pdData = cArray.pdData;
				int iBase = cArray.iBase;
				int iStep = cArray.iStep;
				for(int c = 0;	c < iNrOfCells; c++)
				{
					p3DfFeatureVectors[c][f] = (float)pdData[iBase + c * iStep];
				}
			}

			// create a vector field by passing the pointer to the data represented by p3dv3Data
			ITL_field_regular<VECTOR3> *vectorField = new ITL_field_regular<VECTOR3>(
				&p3DfFeatureVectors[0],
				3,
				pfBlockDimLow,
				pfBlockDimUp,
				piLowPad,
				piHighPad,
				piNeighborhood );

			ITL_globalentropy<VECTOR3> *globalEntropyComputerForVector = new ITL_globalentropy<VECTOR3>( vectorField );

			// convert each vector into bin index
			globalEntropyComputerForVector->computeHistogramBinField("vector", iDefaultNrOfBins);

			// compute the entropy
			globalEntropyComputerForVector->computeGlobalEntropyOfField(iDefaultNrOfBins, false);

			// save the block-wise entropy
			fEntropy = globalEntropyComputerForVector->getGlobalEntropy();

			delete globalEntropyComputerForVector;
		} break;
		}
	#else	// MOD-BY-LEETEN 07/22/2011-TO:
	TBuffer3D<float> p3DfFeatureScalars;

	p3DfFeatureScalars.alloc(
		cBoundBlock.piDimLengths[0],
		cBoundBlock.piDimLengths[1],
		cBoundBlock.piDimLengths[2]
		);

	// MOD-BY-LEETEN 07/31/2011-FROM:
		// _CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0]);
	// TO:
	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0], true);
	// MOD-BY-LEETEN 07/31/2011-END

	ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp,
		piLowPad,
		piHighPad,
		piNeighborhood );

	ITL_globalentropy<SCALAR> *globalEntropyComputerForScalar = new ITL_globalentropy<SCALAR>( scalarField );

	// obtain the default range of the random variable, which is especially useful for orientation
	#if 0 // MOD-BY-LEETEN 07/31/2011-FROM:
		double dDefaultMin, dDefaultMax;
		cRandomVariable._GetDefaultRange(dDefaultMin, dDefaultMax);
		double dHistMin, dHistMax;
		dHistMin = max(cRandomVariable.cRange.dMin, dDefaultMin);
		dHistMax = min(cRandomVariable.cRange.dMax, dDefaultMax);
		globalEntropyComputerForScalar->setHistogramRange(dHistMin, dHistMax);

		// convert each vector into bin index
		globalEntropyComputerForScalar->computeHistogramBinField("scalar", iDefaultNrOfBins);

		// compute the entropy
		globalEntropyComputerForScalar->computeGlobalEntropyOfField(iDefaultNrOfBins, false);
	#else	// MOD-BY-LEETEN 07/31/2011-TO:
	int iNrOfBins = cRandomVariable.IGetNrOfBins();
	globalEntropyComputerForScalar->setHistogramRange(0, iNrOfBins);
	globalEntropyComputerForScalar->computeHistogramBinField("scalar", iNrOfBins);

	// compute the entropy
	globalEntropyComputerForScalar->computeGlobalEntropyOfField(ITL_histogram::DEFAULT_NR_OF_BINS, false);
	#endif	// MOD-BY-LEETEN 07/31/2011-END

	// save the block-wise entropy
	float fEntropy = globalEntropyComputerForScalar->getGlobalEntropy();

	delete globalEntropyComputerForScalar;
	#endif	// MOD-BY-LEETEN 07/22/2011-END

	// save the computed entropy...
	if( NULL != szEntropyLogPathFilename )
	{
		char szCommmand[1024];
		sprintf(szCommmand, "echo Block %d RV %d Entropy %e >> %s", iBoundBlock, iRandomVariable, fEntropy, szEntropyLogPathFilename);
		system(szCommmand);
	}
}

void
ITLRandomField::_ComputeEntorpyFieldInBoundBlock
(
	const int iRandomVariable,
	const int iDim,
	const double pdLocalNeighborhood[],
	const char *szEntropyPathPrefix
)
{
	const int iBoundBlock = IGetBoundBlock();
	const CBlock& cBoundBlock = CGetBoundBlock();

	CRandomVariable& cRandomVariable = CGetRandomVariable(iRandomVariable);
	// DEL-BY-LEETEN 08/06/2011-BEGIN
		// const int iFeatureLength = cRandomVariable.piFeatureVector.USize();
	// DEL-BY-LEETEN 08/06/2011-END

	// compute the #cells
	int iNrOfCells = 1;
	for(int d = 0;	d < CBlock::MAX_DIM; d++)
	{
		iNrOfCells *= cBoundBlock.piDimLengths[d];
	}

	// other setting of the blocks.
	float pfBlockDimLow[CBlock::MAX_DIM];
	float pfBlockDimUp[CBlock::MAX_DIM];
	int	piLowPad[CBlock::MAX_DIM];
	int piHighPad[CBlock::MAX_DIM];
	int piNeighborhood[CBlock::MAX_DIM];
	// DEL-BY-LEETEN 07/31/2011-BEGIN
		// int iDefaultNrOfBins = 360;
	// DEL-BY-LEETEN 07/31/2011-END

	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		pfBlockDimLow[d] = 0.0f;
		pfBlockDimUp[d] = (float)cBoundBlock.piDimLengths[d] - 1;
		piLowPad[d] = 0;
		piHighPad[d] = 0;
		piNeighborhood[d] = (d < iDim)?(int)pdLocalNeighborhood[d]:0;
	}

	#if	0	// MOD-BY-LEETEN 07/22/2011-FROM:
		switch(iFeatureLength)
		{
		case 1:
		{
			TBuffer3D<float> p3DfFeatureScalars;

			p3DfFeatureScalars.alloc(
				cBoundBlock.piDimLengths[0],
				cBoundBlock.piDimLengths[1],
				cBoundBlock.piDimLengths[2]
				);

			for(int f = 0; f < iFeatureLength; f++)
			{
				int iDataComponent = cRandomVariable.piFeatureVector[f];
				CDataComponent& cDataComponent = this->CGetDataComponent(iDataComponent);

				const CArray &cArray = this->CGetArray(iBoundBlock, iDataComponent);
				const double *pdData = cArray.pdData;
				int iBase = cArray.iBase;
				int iStep = cArray.iStep;
				for(int c = 0;	c < iNrOfCells; c++)
				{
					double dValue = (float)pdData[iBase + c * iStep];
					p3DfFeatureScalars[c] = cDataComponent.cRange.DClamp(dValue);
				}
			}

			ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
				&p3DfFeatureScalars[0],
				3,
				pfBlockDimLow,
				pfBlockDimUp,
				piLowPad,
				piHighPad,
				piNeighborhood );

			ITL_localentropy<SCALAR> *localEntropyComputerForScalar = new ITL_localentropy<SCALAR>( scalarField );

			// convert each vector into bin index
			localEntropyComputerForScalar->computeHistogramBinField("scalar", iDefaultNrOfBins);

			// compute the entropy
			localEntropyComputerForScalar->computeEntropyOfField( iDefaultNrOfBins, false);

			if( NULL != szEntropyPathPrefix )
			{
				ITL_field_regular<float>* entropyField = localEntropyComputerForScalar->getEntropyField();
				ITL_ioutil<float>::writeFieldBinarySerial(
						entropyField->getDataFull(),
						szEntropyPathPrefix,
						entropyField->grid->dim,
						3 );
			}

			delete localEntropyComputerForScalar;
		} break;

		case 3:
		{
			TBuffer3D<VECTOR3> p3DfFeatureVectors;

			p3DfFeatureVectors.alloc(
				cBoundBlock.piDimLengths[0],
				cBoundBlock.piDimLengths[1],
				cBoundBlock.piDimLengths[2]
				);

			for(int f = 0; f < iFeatureLength; f++)
			{
				const CArray &cArray = this->CGetArray(iBoundBlock, cRandomVariable.piFeatureVector[f]);
				const double *pdData = cArray.pdData;
				int iBase = cArray.iBase;
				int iStep = cArray.iStep;
				for(int c = 0;	c < iNrOfCells; c++)
				{
					p3DfFeatureVectors[c][f] = (float)pdData[iBase + c * iStep];
				}
			}

			// create a vector field by passing the pointer to the data represented by p3dv3Data
			ITL_field_regular<VECTOR3> *vectorField = new ITL_field_regular<VECTOR3>(
				&p3DfFeatureVectors[0],
				3,
				pfBlockDimLow,
				pfBlockDimUp,
				piLowPad,
				piHighPad,
				piNeighborhood );

			ITL_localentropy<VECTOR3> *localEntropyComputerForVector = new ITL_localentropy<VECTOR3>( vectorField );

			// convert each vector into bin index
			localEntropyComputerForVector->computeHistogramBinField("vector", iDefaultNrOfBins);

			// compute the entropy
			localEntropyComputerForVector->computeEntropyOfField( iDefaultNrOfBins, false);

			if( NULL != szEntropyPathPrefix )
			{
				ITL_field_regular<float>* entropyField = localEntropyComputerForVector->getEntropyField();
				ITL_ioutil<float>::writeFieldBinarySerial(
						entropyField->getDataFull(),
						szEntropyPathPrefix,
						entropyField->grid->dim,
						3 );
			}
			delete localEntropyComputerForVector;
		} break;
		}
	#else	// MOD-BY-LEETEN 07/22/2011-TO:
	TBuffer3D<float> p3DfFeatureScalars;

	p3DfFeatureScalars.alloc(
		cBoundBlock.piDimLengths[0],
		cBoundBlock.piDimLengths[1],
		cBoundBlock.piDimLengths[2]
		);

	// MOD-BY-LEETEN 07/31/2011-FROM:
		// _CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0]);
	// TO:
	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0], true);
	// MOD-BY-LEETEN 07/31/2011-END

	ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp,
		piLowPad,
		piHighPad,
		piNeighborhood );

	ITL_localentropy<SCALAR> *localEntropyComputerForScalar = new ITL_localentropy<SCALAR>( scalarField );

	// obtain the default range of the random variable, which is especially useful for orientation
	#if 0 // MOD-BY-LEETEN 07/31/2011-FROM:
		double dDefaultMin, dDefaultMax;
		cRandomVariable._GetDefaultRange(dDefaultMin, dDefaultMax);
		double dHistMin, dHistMax;
		dHistMin = max(cRandomVariable.cRange.dMin, dDefaultMin);
		dHistMax = min(cRandomVariable.cRange.dMax, dDefaultMax);
		localEntropyComputerForScalar->setHistogramRange(dHistMin, dHistMax);

		// convert each vector into bin index
		localEntropyComputerForScalar->computeHistogramBinField("scalar", iDefaultNrOfBins);

		// compute the entropy
		localEntropyComputerForScalar->computeEntropyOfField( iDefaultNrOfBins, false);
	#else	// MOD-BY-LEETEN 07/31/2011-TO:
	int iNrOfBins = cRandomVariable.IGetNrOfBins();
	localEntropyComputerForScalar->setHistogramRange(0, iNrOfBins);

	// convert each vector into bin index
	localEntropyComputerForScalar->computeHistogramBinField("scalar", iNrOfBins);

	// compute the entropy
	localEntropyComputerForScalar->computeEntropyOfField( ITL_histogram::DEFAULT_NR_OF_BINS, false);
	#endif	// MOD-BY-LEETEN 07/31/2011-END

	if( NULL != szEntropyPathPrefix )
	{
		ITL_field_regular<float>* entropyField = localEntropyComputerForScalar->getEntropyField();
		ITL_ioutil<float>::writeFieldBinarySerial(
				entropyField->getDataFull(),
				szEntropyPathPrefix,
				entropyField->grid->dim,
				3 );
	}

	delete localEntropyComputerForScalar;
	#endif	// MOD-BY-LEETEN 07/22/2011-END
}

// ADD-BY-LEETEN 07/31/2011-BEGIN
//////////////////////////////////////////////////////////////////////////
void
ITLRandomField::_ComputeJointEntorpyInBoundBlock
(
	int iRandomVariable1,
	int iRandomVariable2,
	const char *szEntropyLogPathFilename
)
{
	const int iBoundBlock = IGetBoundBlock();
	const CBlock& cBoundBlock = CGetBoundBlock();

	CRandomVariable& cRandomVariable1 = CGetRandomVariable(iRandomVariable1);
	CRandomVariable& cRandomVariable2 = CGetRandomVariable(iRandomVariable2);

	// compute the #cells
	int iNrOfCells = 1;
	for(int d = 0;	d < CBlock::MAX_DIM; d++)
	{
		iNrOfCells *= cBoundBlock.piDimLengths[d];
	}

	// other setting of the blocks.
	float pfBlockDimLow[CBlock::MAX_DIM];
	float pfBlockDimUp[CBlock::MAX_DIM];
	int	piLowPad[CBlock::MAX_DIM];
	int piHighPad[CBlock::MAX_DIM];
	int piNeighborhood[CBlock::MAX_DIM];

	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		pfBlockDimLow[d] = 0.0f;
		pfBlockDimUp[d] = (float)cBoundBlock.piDimLengths[d] - 1;
		piLowPad[d] = 0;
		piHighPad[d] = 0;
		piNeighborhood[d] = 0;
	}

	TBuffer3D<float> p3DfFeatureScalars1;
	p3DfFeatureScalars1.alloc(
		cBoundBlock.piDimLengths[0],
		cBoundBlock.piDimLengths[1],
		cBoundBlock.piDimLengths[2]
		);
	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable1, &p3DfFeatureScalars1[0], true);

	ITL_field_regular<SCALAR> *scalarField1 = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars1[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp,
		piLowPad,
		piHighPad,
		piNeighborhood );

	TBuffer3D<float> p3DfFeatureScalars2;
	p3DfFeatureScalars2.alloc(
		cBoundBlock.piDimLengths[0],
		cBoundBlock.piDimLengths[1],
		cBoundBlock.piDimLengths[2]
		);
	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable2, &p3DfFeatureScalars2[0], true);

	ITL_field_regular<SCALAR> *scalarField2 = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars2[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp,
		piLowPad,
		piHighPad,
		piNeighborhood );

	// Initialize class that can compute global joint entropy
	ITL_globaljointentropy<SCALAR>* globalJointEntropyComputer = new ITL_globaljointentropy<SCALAR>( scalarField1, scalarField2 );

	int iNrOfBins1 = cRandomVariable1.IGetNrOfBins();
	int iNrOfBins2 = cRandomVariable2.IGetNrOfBins();
	globalJointEntropyComputer->setJointHistogramRange( 0, iNrOfBins1, 0, iNrOfBins2 );	
	globalJointEntropyComputer->computeJointHistogramBinField( "scalar",  ITL_histogram::DEFAULT_NR_OF_BINS);

	// Global Joint entropy Computation
	globalJointEntropyComputer->computeGlobalJointEntropyOfField( ITL_histogram::DEFAULT_NR_OF_BINS, false );
	
	// save the block-wise entropy
	float fJointEntropy = globalJointEntropyComputer->getGlobalJointEntropy();

	// save the computed entropy...
	if( NULL != szEntropyLogPathFilename )
	{
		char szCommmand[1024];
		sprintf(szCommmand, "echo Block %d RV %d %d JEntropy %e >> %s", iBoundBlock, iRandomVariable1, iRandomVariable2, fJointEntropy, szEntropyLogPathFilename);
		system(szCommmand);
	}

	delete globalJointEntropyComputer;
}
// ADD-BY-LEETEN 07/31/2011-END

//////////////////////////////////////////////////////////////////////////
ITLRandomField::ITLRandomField() {
	// TODO Auto-generated constructor stub
	this->iBoundBlock = -1;
	this->iBoundDataComponent = -1;
	this->iBoundRandomVariable = -1;
	// ADD-BY-LEETEN 08/06/2011-BEGIN
	iNcId = 0;
	iRank = -1;
	memset(piNcDimIds, 0, sizeof(piNcDimIds));
	memset(piNcDimVarIds, 0, sizeof(piNcDimVarIds));
	iNcTimeVarId = 0;
	viTimeStamps.push_back(GLOBAL_TIME_STAMP);
	// ADD-BY-LEETEN 08/06/2011-END
	iNrOfGlobalBlocks = 0;	// ADD-BY-LEETEN 08/12/2011
}

ITLRandomField::~ITLRandomField() {
	// TODO Auto-generated destructor stub
  // ADD-BY-LEETEN 07/23/2011-BEGIN
  for(vector<CRandomVariable*>::iterator ivcRandomVariables = vcRandomVariables.begin();
      ivcRandomVariables != vcRandomVariables.end();
      ivcRandomVariables++)
	if( *ivcRandomVariables )
	{
		delete *ivcRandomVariables;
		*ivcRandomVariables = NULL;
	}
  // ADD-BY-LEETEN 07/23/2011-BEGIN
}
