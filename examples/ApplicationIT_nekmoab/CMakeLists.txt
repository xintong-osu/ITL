# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (NekMoab) 

# find the mpi package
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)

# use mpicxx to replace 
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# specify the path to the headers
include_directories(
		. 
		./include
		${ITL_DIR}/include
	)
include_directories( 
		${DIY_INC}
	)
include_directories( 
		${BIL_INC}
	)
include_directories( 
		${MOAB_INC}
	)
include_directories( 
		${VTK_INC}
	)

# grab all the source codes
AUX_SOURCE_DIRECTORY ("src" 	SRC_FILES)

# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}_r"
)

install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)
#
# $Log$
#

