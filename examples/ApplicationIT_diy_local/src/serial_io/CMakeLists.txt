cmake_minimum_required( VERSION 2.8 )
if( WIN32 )
	cmake_policy( SET CMP0015 NEW )
endif()

project( serialIO )

# Add third party dependencies
find_package(MPI REQUIRED)

# Use MPI compiler and linker
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Use mpicxx to replace 
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

################################################
# Add third party dependencies
################################################
# Set DIY path
set(DIY_DIR "../.." CACHE PATH "The path to DIY")

# Set zlib path
set(ZLIB_DIR "../.." CACHE PATH "The path to ZLIB")

################################################
# Specify the path to the headers
################################################
# Include various third party headers
# MPI
include_directories( "${MPI_INCLUDE_PATH}" )
#####
# Our libraries
# DIY
include_directories( "${DIY_DIR}/include" )
#
include_directories( "${PROJECT_SOURCE_DIR}" )

################################################
# Specify path to the various libraries to link
################################################
list(APPEND LIB_DIRS "${MPI_DIR}/lib" )
#####
list(APPEND LIB_DIRS "${ZLIB_DIR}/lib" )
#####
# Our libraries
# DIY
list(APPEND LIB_DIRS "${DIY_DIR}/lib" )

################################################
# Link libraries
################################################
link_directories( ${LIB_DIRS} )

################################################
# Build the project as an executable 
# Grab all the source codes
################################################ 
add_executable( serialIO ${PROJECT_SOURCE_DIR}/Main_reader.cpp )

################################################
# link the libraries to the executable
################################################
#####
if( WIN32 )
    target_link_libraries( ${PROJECT_NAME} mpi )
else()
    target_link_libraries( ${PROJECT_NAME} mpich )
endif()
#####
target_link_libraries( ${PROJECT_NAME} z )
#####
target_link_libraries( ${PROJECT_NAME} diy )
#####

# Installation Instructions
# bin is interpreted as relative to CMAKE_INSTALL_PREFIX
install( TARGETS ${PROJECT_NAME} DESTINATION bin )


