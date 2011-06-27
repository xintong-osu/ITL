# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (ApplicationIT_regularvector) 

# set ITL path
set(ITL_DIR "../.." CACHE PATH "The path to ITL")

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
		${ITL_DIR}/include/datastore
		${ITL_DIR}/include/field
		${ITL_DIR}/include/grid
		${ITL_DIR}/include/metrics
		${ITL_DIR}/include/util
	)

link_directories(
	${ITL_DIR}/lib
)

# grab all the source codes
AUX_SOURCE_DIRECTORY ("src" 	SRC_FILES)

# build the project as a library
add_executable(${PROJECT_NAME} ${SRC_FILES})

# link the libraries
target_link_libraries(${PROJECT_NAME} ITLib) 

#
# $Log$
#

