/*
 * ITL_random_field.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: leeten
 */

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
	// MOD-BY-LEETEN 09/01/2011-FROM:
		// ASSERT_OR_LOG(MPI_SUCCESS == MPI_Allreduce(MPI_IN_PLACE, &iNrOfGlobalBlocks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD), "");
	// TO:
	ASSERT_OR_LOG(MPI_SUCCESS == MPI_Allreduce(&iMaxGlobalId, &iNrOfGlobalBlocks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD), "");
	// MOD-BY-LEETEN 09/01/2011-END
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

#if	0	// MOD-BY-LEETEN 08/29/2011-FROM:
	float fDomainMax;
	MPI_Reduce(&fLocalMax, &fDomainMax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	float fDomainMin;
	MPI_Reduce(&fLocalMin, &fDomainMin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

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
	// ADD-BY-LEETEN 09/01/2011-BEGIN
	if( iNcId >= 0 )
	{
		// write the time stamp
		#ifndef	WITH_PNETCDF		// ADD-BY-LEETEN 08/12/2011
	        size_t uStart = viTimeStamps.size();
		size_t uCount = 1;
		ASSERT_NETCDF(nc_put_vara_int(
				iNcId,
				iNcTimeVarId,
				&uStart,
				&uCount,
				&iTimeStamp));

		// ADD-BY-LEETEN 08/12/2011-BEGIN
		#else	// #ifndef	WITH_PNETCDF
		MPI_Offset uStart = viTimeStamps.size();
		MPI_Offset uCount = 1;

		ASSERT_NETCDF(ncmpi_begin_indep_data(iNcId));
		if( 0 == iRank )
			ASSERT_NETCDF(ncmpi_put_vara_int(
					iNcId,
					iNcTimeVarId,
					&uStart,
					&uCount,
					&iTimeStamp));
		ASSERT_NETCDF(ncmpi_end_indep_data(iNcId));
		#endif	// #ifndef	WITH_PNETCDF
		// ADD-BY-LEETEN 08/12/2011-END
	}
	// ADD-BY-LEETEN 09/01/2011-END
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
	vcRandomVariables.push_back(new CRandomVariable());
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
	return *vcRandomVariables[iRandomVariable];
}

const ITLRandomField::CRandomVariable&
ITLRandomField::CGetRandomVariable
(
	const int iRandomVariable
) const
{
	return *vcRandomVariables[iRandomVariable];
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
	const int iFeatureMapping
)
{
	CGetBoundRandomVariable()._Set(iFeatureLength, piFeatureVector, iFeatureMapping);
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
		dSample = cRandomVariable.cRange.DClamp(dSample);
		if( bIsMappingToBins )
		  {
		    dSample = cRandomVariable.DMapToBin(dSample);
		  }
		pfRandomSamples[c] = (float)dSample;
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
	this->pcDataComponentInfo.New(iNrOfDataComponents);

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
	// MOD-BY-LEETEN 09/01/2011-FROM:
		// if( iNcId > 0 )
	// TO:
	if( iNcId >= 0 )
	// MOD-BY-LEETEN 09/01/2011-END
	{
		// write the time stamp
		TBuffer<int> piTemp;
		piTemp.alloc(this->IGetNrOfTimeStamps());
		for(int t = 0; t < (int)piTemp.USize(); t++)
			piTemp[t] = this->viTimeStamps[t];

		#ifndef	WITH_PNETCDF		// ADD-BY-LEETEN 08/12/2011
		#if 0 // DEL-BY-LEETEN 09/01/2011-BEGIN
			// since the time step wil lbe written earlier, this part can be removed
			size_t uStart = 0;
			size_t uCount = piTemp.USize();
			ASSERT_NETCDF(nc_put_vara_int(
					iNcId,
					iNcTimeVarId,
					&uStart,
					&uCount,
					&piTemp[0]));
		#endif	// DEL-BY-LEETEN 09/01/2011-END
        /* Close the file. */
	    ASSERT_NETCDF(nc_close(iNcId));

		// ADD-BY-LEETEN 08/12/2011-BEGIN
		#else	// #ifndef	WITH_PNETCDF

	    #if 0 // DEL-BY-LEETEN 09/01/2011-BEGIN
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
		#endif	// DEL-BY-LEETEN 09/01/2011-END

        /* Close the file. */
	    ASSERT_NETCDF(ncmpi_close(iNcId));
		#endif	// #ifndef	WITH_PNETCDF
		// ADD-BY-LEETEN 08/12/2011-END

	    // MOD-BY-LEETEN 09/01/2011-FROM:
	    	// iNcId = 0;
		// TO:
	    iNcId = -1;
	    // MOD-BY-LEETEN 09/01/2011-END
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
		#ifndef WITH_PNETCDF		// ADD-BY-LEETEN 09/01/2011
		size_t puStart[3];
		size_t puCount[3];

		// ADD-BY-LEETEN 09/01/2011-BEGIN
		#else	// #ifndef WITH_PNETCDF
		MPI_Offset puStart[3];
		MPI_Offset puCount[3];
		MPI_Offset puStride[3];
		for(int i = 0; i < sizeof(puStride)/sizeof(puStride[0]); i++)
			puStride[i] = 1;
		#endif	// #ifndef WITH_PNETCDF
		// ADD-BY-LEETEN 09/01/2011-END

		// time
		puStart[0] = this->IGetNrOfTimeStamps() - 1;
		puCount[0] = 1;

		// block
		#ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
		puStart[1] = (size_t)iBlockId;
		// ADD-BY-LEETEN 09/01/2011-BEGIN
		#else	// #ifndef WITH_PNETCDF
		puStart[1] = (size_t)CGetBlock(iBlockId).iGlobalId;
		#endif	// #ifndef WITH_PNETCDF
		// ADD-BY-LEETEN 09/01/2011-END
		puCount[1] = 1;

		// spatial
		puStart[2] = 0;
		puCount[2] = (size_t)this->CGetBlock(iBlockId).piDimLengths[d];

		this->CGetBlock(iBlockId)._DumpDimGeometry2Nc(
						d,
						iNcId,
						piNcDimVarIds[d],
						#ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
						puStart,
						puCount);
						// ADD-BY-LEETEN 09/01/2011-BEGIN
						#else	// #ifndef WITH_PNETCDF
						puStart,
						puCount,
						puStride);
						#endif	// #ifndef WITH_PNETCDF
						// ADD-BY-LEETEN 09/01/2011-END
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

		#ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
      size_t puStart[6];
      size_t puCount[6];
      // ADD-BY-LEETEN 09/01/2011-BEGIN
      #else	// #ifndef WITH_PNETCDF
      MPI_Offset puStart[6];
      MPI_Offset puCount[6];
      MPI_Offset puStride[6];
      for(int d = 0; d < sizeof(puStride)/sizeof(puStride[0]); d ++)
    	  puStride[d] = 1;
		#endif	// #ifndef WITH_PNETCDF
      // ADD-BY-LEETEN 09/01/2011-END

      // time
      puStart[0] = this->IGetNrOfTimeStamps() - 1;
      puCount[0] = 1;

      // block
		#ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
      puStart[1] = (size_t)b;
      // ADD-BY-LEETEN 09/01/2011-BEGIN
      #else	// #ifndef WITH_PNETCDF
      puStart[1] = (size_t)CGetBlock(b).iGlobalId;
	  #endif	// #ifndef WITH_PNETCDF
      // ADD-BY-LEETEN 09/01/2011-END
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
	  	// MOD-BY-LEETEN 09/01/2011-FROM:
			// LOG_ERROR(fprintf(stderr, "PNetCDF is not fully supported yet."))
		// TO:
		ASSERT_NETCDF(ncmpi_put_vars_double_all(
					   iNcId,
					   cDataComponent.iVarId,
					   puStart,
					   puCount,
					   puStride,
					   &pdTemp[0]));
		// MOD-BY-LEETEN 09/01/2011-END
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

	  #ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
      size_t puStart[6];
      size_t puCount[6];

      // ADD-BY-LEETEN 09/01/2011-BEGIN
      #else	// #ifndef WITH_PNETCDF
      MPI_Offset puStart[6];
      MPI_Offset puCount[6];
      MPI_Offset puStride[6];
      for(int d = 0; d < sizeof(puStride)/sizeof(puStride[0]); d ++)
    	  puStride[d] = 1;
	  #endif	// #ifndef WITH_PNETCDF
      // ADD-BY-LEETEN 09/01/2011-END

      // time
      puStart[0] = this->IGetNrOfTimeStamps() - 1;
      puCount[0] = 1;

      // block
	  #ifndef WITH_PNETCDF	// ADD-BY-LEETEN 09/01/2011
      puStart[1] = (size_t)b;
      // ADD-BY-LEETEN 09/01/2011-BEGIN
      #else	// #ifndef WITH_PNETCDF
      puStart[1] = (size_t)CGetBlock(b).iGlobalId;
	  #endif	// #ifndef WITH_PNETCDF
      // ADD-BY-LEETEN 09/01/2011-END
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
      	// MOD-BY-LEETEN 09/01/2011-FROM:
			// LOG_ERROR(fprintf(stderr, "PNetCDF is not fully supported yet."))
		// TO:
		ASSERT_NETCDF(ncmpi_put_vars_float_all(
					  iNcId,
					  CGetRandomVariable(iRandomVariable).iVarId,
					  puStart,
					  puCount,
					  puStride,
					  &pfTemp[0]));
		// MOD-BY-LEETEN 09/01/2011-END
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

	for(int d = 0; d < CBlock::MAX_DIM; d++)
	{
		pfBlockDimLow[d] = 0.0f;
		pfBlockDimUp[d] = (float)cBoundBlock.piDimLengths[d] - 1;
		piLowPad[d] = 0;
		piHighPad[d] = 0;
		piNeighborhood[d] = 0;
	}

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

	#if	0	// MOD-BY-LEETEN 12/02/2011-FROM:
		ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
			&p3DfFeatureScalars[0],
			3,
			pfBlockDimLow,
			pfBlockDimUp,
			piLowPad,
			piHighPad,
			piNeighborhood );
	#else	// MOD-BY-LEETEN 12/02/2011-TO:
	ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp);
	#endif	// MOD-BY-LEETEN 12/02/2011-END

	// MOD-BY-LEETEN 12/02/2011-FROM:
		// ITL_globalentropy<SCALAR> *globalEntropyComputerForScalar = new ITL_globalentropy<SCALAR>( scalarField );
	// TO:
	ITL_globalentropy<SCALAR> *globalEntropyComputerForScalar = new ITL_globalentropy<SCALAR>( scalarField, pcHistogram);
	// MOD-BY-LEETEN 12/02/2011-END

	// obtain the default range of the random variable, which is especially useful for orientation
	int iNrOfBins = cRandomVariable.IGetNrOfBins();
	globalEntropyComputerForScalar->setHistogramRange(0, iNrOfBins);
	globalEntropyComputerForScalar->computeHistogramBinField("scalar", iNrOfBins);

	// compute the entropy
	globalEntropyComputerForScalar->computeGlobalEntropyOfField(ITL_histogram::DEFAULT_NR_OF_BINS, false);

	// save the block-wise entropy
	float fEntropy = globalEntropyComputerForScalar->getGlobalEntropy();

	delete globalEntropyComputerForScalar;

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
		piNeighborhood[d] = (d < iDim)?(int)pdLocalNeighborhood[d]:0;
	}

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

	#if	0	// MOD-BY-LEETEN 12/02/2011-FROM:
		ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
			&p3DfFeatureScalars[0],
			3,
			pfBlockDimLow,
			pfBlockDimUp,
			piLowPad,
			piHighPad,
			piNeighborhood );
	#else	// MOD-BY-LEETEN 12/02/2011-TO:
	ITL_field_regular<SCALAR> *scalarField = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp);
	#endif	// MOD-BY-LEETEN 12/02/2011-END

	// MOD-BY-LEETEN 12/02/2011-FROM:
		// ITL_localentropy<SCALAR> *localEntropyComputerForScalar = new ITL_localentropy<SCALAR>( scalarField );
	// TO:
	ITL_localentropy<SCALAR> *localEntropyComputerForScalar = new ITL_localentropy<SCALAR>( scalarField, pcHistogram );
	// MOD-BY-LEETEN 12/02/2011-END

	// obtain the default range of the random variable, which is especially useful for orientation
	int iNrOfBins = cRandomVariable.IGetNrOfBins();
	localEntropyComputerForScalar->setHistogramRange(0, iNrOfBins);

	// convert each vector into bin index
	localEntropyComputerForScalar->computeHistogramBinField("scalar", iNrOfBins);

	// compute the entropy
	localEntropyComputerForScalar->computeEntropyOfField( ITL_histogram::DEFAULT_NR_OF_BINS, false);

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

	#if	0	// MOD-BY-LEETEN 12/02/2011-FROM:
		ITL_field_regular<SCALAR> *scalarField1 = new ITL_field_regular<SCALAR>(
			&p3DfFeatureScalars1[0],
			3,
			pfBlockDimLow,
			pfBlockDimUp,
			piLowPad,
			piHighPad,
			piNeighborhood );
	#else	// MOD-BY-LEETEN 12/02/2011-TO:
	ITL_field_regular<SCALAR> *scalarField1 = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars1[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp);
	#endif	// MOD-BY-LEETEN 12/02/2011-END

	TBuffer3D<float> p3DfFeatureScalars2;
	p3DfFeatureScalars2.alloc(
		cBoundBlock.piDimLengths[0],
		cBoundBlock.piDimLengths[1],
		cBoundBlock.piDimLengths[2]
		);
	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable2, &p3DfFeatureScalars2[0], true);

	#if	0	// MOD-BY-LEETEN 12/02/2011-FROM:
		ITL_field_regular<SCALAR> *scalarField2 = new ITL_field_regular<SCALAR>(
			&p3DfFeatureScalars2[0],
			3,
			pfBlockDimLow,
			pfBlockDimUp,
			piLowPad,
			piHighPad,
			piNeighborhood );
	#else	// MOD-BY-LEETEN 12/02/2011-TO:
	ITL_field_regular<SCALAR> *scalarField2 = new ITL_field_regular<SCALAR>(
		&p3DfFeatureScalars2[0],
		3,
		pfBlockDimLow,
		pfBlockDimUp);
	#endif	// MOD-BY-LEETEN 12/02/2011-END

	// Initialize class that can compute global joint entropy
	// MOD-BY-LEETEN 12/02/2011-FROM:
		// ITL_globaljointentropy<SCALAR>* globalJointEntropyComputer = new ITL_globaljointentropy<SCALAR>( scalarField1, scalarField2 );
	// TO:
	ITL_globaljointentropy<SCALAR>* globalJointEntropyComputer = new ITL_globaljointentropy<SCALAR>( scalarField1, scalarField2, pcHistogram );
	// MOD-BY-LEETEN 12/02/2011-END

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
	// MOD-BY-LEETEN 09/01/2011-FROM:
		// iNcId = 0;
	// TO:
	iNcId = -1;
	// MOD-BY-LEETEN 09/01/2011-END
	iRank = -1;
	memset(piNcDimIds, 0, sizeof(piNcDimIds));
	memset(piNcDimVarIds, 0, sizeof(piNcDimVarIds));
	iNcTimeVarId = 0;
	viTimeStamps.push_back(GLOBAL_TIME_STAMP);
	// ADD-BY-LEETEN 08/06/2011-END
	iNrOfGlobalBlocks = 0;	// ADD-BY-LEETEN 08/12/2011

	// ADD-BY-LEETEN 12/02/2011-BEGIN
	pcHistogram = NULL;

	// Initialize histogram
	pcHistogram = new ITL_histogram( "!" );
	// ADD-BY-LEETEN 12/02/2011-END
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

	// ADD-BY-LEETEN 12/02/2011-BEGIN
  	  if( pcHistogram )
  	  {
  		delete pcHistogram;
  		pcHistogram = NULL;
  	  }
	// ADD-BY-LEETEN 12/02/2011-END
}
