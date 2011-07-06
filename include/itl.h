/*---------------------------------------------------------------------------
 *
 * itl C, C++ interface
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * Copyright Notice
 * + 2010 University of Chicago
 *
--------------------------------------------------------------------------*/

#ifndef _ITL
#define _ITL

#if	0	// MOD-BY-LEETEN 07/05/2011-FROM:
	#include "mpi.h"

	/*-------------------------------------------------------------------------*/

	#ifdef __cplusplus
	extern "C"
	#endif
	void ITL_begin();
	#ifdef __cplusplus
	extern "C"
	#endif
	void ITL_get_data();

	/*-------------------------------------------------------------------------*/
#else	// MOD-BY-LEETEN 07/05/2011-TO:

//--------------------------------------------------------------------------
// functions

//! The C and C++ API to initialize ITL
/*!
 *
*/
void
ITL_begin();

//! The Fortran API to initialize ITL
/*!
 * \sa ITL_begin
*/
extern "C"
void
itl_begin_();

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to free ITL
/*!
 *
*/
void
ITL_end();

//! The Fortran API to free ITL
/*!
 * \sa ITL_end
*/
extern "C"
void itl_end_();

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify #blocks
/*!
\param iNrOfBlocks		#blocks
*/
void
ITL_nblocks
(
	const int iNrOfBlocks
);

//! The Fortran API to specify #blocks and length of each element
/*!
 * \sa ITL_nblocks
*/
extern "C"
void
itl_nblocks_
(
	int *piNrOfBlocks
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify length of feature vector on each element
/*!
\param iFeatureLength 	length of each elements. Currently only 1, 2, and 3 are supported.
*/
void
ITL_feature_length
(
	const int iFeatureLength
);

//! The Fortran API to specify length of feature vector on each element
/*!
 * \sa ITL_feature_length
*/
extern "C"
void
itl_feature_length_
(
	int *piFeatureLength
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to bind the block.
/*!
\param iBlockId Id of the block (0-based)
*/
void
ITL_bind_block
(
	const int iBlockId
);

//! The Fortran API to bind the block.
/*!
\param piBlockId pointer to the block id (1-based)
\sa ITL_bind_block
*/
extern "C"
void
itl_bind_block_
(
	int *piBlockId
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the block order
/*!
\param iOrder	Order of the data block.
				0: row-major order.
				1: column-major order
				2: no order (treat it as a 1D array)
				NOTE: currenly only row-major order is supported
*/
void
ITL_block_order
(
	const int iOrder
);

//! The Fortran API to specify the block order
/*!
\sa   ITL_block_order
*/
extern "C"
void
itl_block_order_
(
	int *piOrder
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the block size
/*!
\param iBlockDim 	Dimension of the block
\param piBlockDimLengths
					The length of each dimension. The size should be
					equal to iBlockDim. If length is large than MAX_BLOCK_DIM,
					those elements exceed MAX_BLOCK_DIM will be ignored.
\sa   MAX_BLOCK_DIM
*/
void
ITL_block_size
(
	const int iBlockDim,
	const int piBlockDimLengths[]
);

//! The Fortran API to specify the block size
/*!
\sa   ITL_block_size
*/
extern "C"
void
itl_block_size_
(
	int *piBlockDim,
	int *piBlockDimLengths
);

//! The Fortran API to specify the size for 2D block
/*!
\sa   ITL_block_size
*/
extern "C"
void
itl_block_size2_
(
	int *piBlockXLength,
	int *piBlockYLength
);

//! The Fortran API to specify the size for 3D block
/*!
\sa   ITL_block_size
*/
extern "C"
void
itl_block_size3_
(
	int *piBlockXLength,
	int *piBlockYLength,
	int *piBlockZLength
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the extent within the block.
/*!
\param iBlockDim 	Dimension of the block
\param piBlockDimLow	Lower bound of each dimension
\param piBlockDimUp		Upper bound of each dimension
					The extent of each dimension. The unit is in cells.
					The size should be equal to iBlockDim. If length is
					large than MAX_BLOCK_DIM, those elements exceed MAX_BLOCK_DIM will be ignored.

\sa   MAX_BLOCK_DIM

The C and C++ API to specify the extend within the block. The entropy will
be only computed within that region.
*/
void
ITL_block_extent_in_cells
(
	const int iBlockDim,
	const int piBlockDimLow[],
	const int piBlockDimUp[]
);

//! The Fortran API to specify the extent within the block.
/*!
\sa   ITL_block_extent_in_cells
*/
extern "C"
void
itl_block_extent_in_cells_
(
	int *piBlockDim,
	int *piBlockDimLow,
	int *piBlockDimUp
);

//! The Fortran API to specify the 2D extent within the block.
/*!
\sa   ITL_block_extent_in_cells
*/
extern "C"
void
itl_block_extent_in_cells2_
(
	int *piBlockXLow,
	int *piBlockYLow,
	int *piBlockXUp,
	int *piBlockYUp
);

//! The Fortran API to specify the 3D extent within the block.
/*!
\sa   ITL_block_extent_in_cells
*/
extern "C"
void
itl_block_extent_in_cells3_
(
	int *piBlockXLow,
	int *piBlockYLow,
	int *piBlockZLow,
	int *piBlockXUp,
	int *piBlockYUp,
	int *piBlockZUp
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the neighborhood size
/*!
\param iBlockDim 	Dimension of the block
\param pdLocalNeighborhood
					The length of each dimension. The size should be
					equal to iBlockDim. If length is large than MAX_BLOCK_DIM,
					those elements exceed MAX_BLOCK_DIM will be ignored.
\sa   MAX_BLOCK_DIM
*/
void
ITL_local_neighborhood_in_cells
(
	const int iBlockDim,
	const double pdLocalNeighborhood[]
);

//! The Fortran API to specify the neighborhood size
/*!
\sa   ITL_local_neighborhood_in_cells
*/
extern "C"
void
itl_local_neighborhood_in_cells_
(
	int *piBlockDim,
	double *pdNeighborhood
);

//! The Fortran API to specify the 2D neighborhood size
/*!
\sa   ITL_local_neighborhood_in_cells
*/
extern "C"
void
itl_local_neighborhood2_
(
	double *pdXNeighborhood,
	double *pdYNeighborhood
);

//! The Fortran API to specify the 3D neighborhood size
/*!
\sa   ITL_local_neighborhood_in_cells
*/
extern "C"
void
itl_local_neighborhood3_
(
	double *pdXNeighborhood,
	double *pdYNeighborhood,
	double *pdZNeighborhood
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the pointer, start, and step to access the data
/*!
\param	iFeatureIndex	Index (0-based) of the feature vector
\param	pdData			The data pool
\param	iBase 			Index (0-based) to the 1st element in the pool
\param	iStep			Index difference between consecutive elements in the pool
*/
void
ITL_feature_data
(
	const int iFeatureIndex,
	double pdData[],
	const int iBase,
	const int iStep
);

//! The Fortran ++ API to specify the pointer, start, and step to access the data
/*!
\param	piFeatureIndex	Pointer to the index (1-based) of the feature vector
\param	pdData			The data pool
\param	iBase 			Pointer to the index (1-based) to the 1st element in the pool
\param	iStep			Pointer to the index difference between consecutive elements in the pool

\sa		ITL_feature_data
*/
extern "C"
void
itl_feature_data_
(
	int *piFeatureIndex,
	double *pdData,
	int *piBase,
	int *piStep
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the vector orientation is used as the random variable
/*!
*/
void
ITL_is_using_vector_orientation
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_using_vector_orientation
*/
extern "C"
void
itl_is_using_vector_orientation_
(
	int *piIsEnabled
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the entropy field is computed
/*!
*/
void
ITL_is_computing_local_entropy
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_computing_local_entropy
*/
extern "C"
void
itl_is_computing_local_entropy_
(
	int *piIsComputingLocalEntropy
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to trigger the computation of metrics
/*!
*/
void
ITL_execute
(
);

//! The Fortran API to trigger the computation of metrics
/*!
 * \sa ITL_execute
*/
extern "C"
void
itl_execute_
(
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the geoometry is dumped
/*!
*/
void
ITL_is_dumping_geometry
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_dumping_geometry
*/
extern "C"
void
itl_is_dumping_geometry_
(
	int *piIsEnabled
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the random variable is dumped
/*!
*/
void
ITL_is_dumping_samples
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_dumping_samples
*/
extern "C"
void
itl_is_dumping_samples_
(
	int *piIsEnabled
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the entropy is dumped
/*!
*/
void
ITL_is_dumping_entropy
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_dumping_entropy
*/
extern "C"
void
itl_is_dumping_entropy_
(
	int *piIsEnabled
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify whether the computation is applied to the entire domain.
/*!
 * If true, when call ITL_Execute, it will check where all blocks has been specified.
 * If not, an error message will be thrown.
 * Otherwise, ITL_Execute() will only compute the entropy within the current bound block.
*/
void
ITL_is_executing_all_blocks
(
	bool bIsEnabled
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \sa ITL_is_executing_all_blocks
*/
extern "C"
void
itl_is_executing_all_blocks_
(
	int *piIsEnabled
);

/////////////////////////////////////////////////////////////////////
//! The C and C++ API to specify the coordinates along one dim. for regular grid
/*!
 * \param iDimId	ID (0-based) of the block dimension
 * \param pdCoord	the pool of the coordinates along the specified block dim
 * \param iBase		the 1st element (0-based) in the pool
 * \param iStep		the difference between the consecutive elements
*/
void
ITL_geom_rect_dim_coord
(
	const int iDimId,
	double *pdCoord,
	const int iBase,
	const int iStep
);

//! The Fortran API to specify whether the vector orientation is used as the random variable
/*!
 * \param iDimId	ID (1-based) of the block dimension
 * \param pdCoord	the pool of the coordinates along the specified block dim
 * \param iBase		the 1st element (1-based) in the pool
 * \param iStep		the difference between the consecutive elements
 *
 * \sa	ITL_geom_rect_dim_coord
*/
extern "C"
void
itl_geom_rect_dim_coord_
(
	int *piDimId,
	double *pdCoord,
	int *piBase,
	int *piStep
);

#endif	// MOD-BY-LEETEN 07/05/2011-END
#endif	// #ifndef _ITL

/*
 *
 * $Log$
 *
 */
