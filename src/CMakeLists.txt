# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITLib) 

set(ITL_DIR "/homes/leeten/workspace/itl" CACHE PATH "The path to ITL")
set(MYLIB_DIR "/homes/leeten/mylib" CACHE PATH "The path to MYLIB")

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
		../include
		${ITL_DIR}/include
		${MYLIB_DIR}/include
	)

# grab all the source codes
AUX_SOURCE_DIRECTORY ("." 	SRC_FILES)

# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES})

# output the library to the lib/ folder
set(LIBRARY_OUTPUT_PATH "../lib")

#
# $Log$
#

