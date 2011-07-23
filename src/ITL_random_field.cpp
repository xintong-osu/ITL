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
	#include "ITL_field_regular.h"

	#include "ITL_random_field.h"

// ADD-BY-LEETEN 07/22/2011-BEGIN
void
ITLRandomField::_UseDomainRange
(
	const int iRandomVariable
)
{
	CRandomVariable& cRandomVariable = CGetRandomVariable(iRandomVariable);
	const int iFeatureLength = cRandomVariable.piFeatureVector.USize();

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

	int iRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &iRank);

	float fDomainMax;
	MPI_Reduce(&fLocalMax, &fDomainMax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	float fDomainMin;
	MPI_Reduce(&fLocalMin, &fDomainMin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

	if( 0 == iRank )
	{
		LOG_VAR_TO_ERROR(fDomainMin);
		LOG_VAR_TO_ERROR(fDomainMax);
	}

	MPI_Bcast(&fDomainMin, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&fDomainMax, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	cRandomVariable.cRange._Set((double)fDomainMin, (double)fDomainMax);
}
// ADD-BY-LEETEN 07/22/2011-END

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
	float pfRandomSamples[]
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
		pfRandomSamples[c] = (float)cRandomVariable.cRange.DClamp(dSample);
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
}

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
	const int iFeatureLength = cRandomVariable.piFeatureVector.USize();

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
	int iDefaultNrOfBins = 360;

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

	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0]);

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
	const int iFeatureLength = cRandomVariable.piFeatureVector.USize();

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
	int iDefaultNrOfBins = 360;

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

	_CollectRandomSamplesInBlock(iBoundBlock, iRandomVariable, &p3DfFeatureScalars[0]);

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

//////////////////////////////////////////////////////////////////////////
ITLRandomField::ITLRandomField() {
	// TODO Auto-generated constructor stub
	this->iBoundBlock = -1;
	this->iBoundDataComponent = -1;
	this->iBoundRandomVariable = -1;
}

ITLRandomField::~ITLRandomField() {
	// TODO Auto-generated destructor stub
}
