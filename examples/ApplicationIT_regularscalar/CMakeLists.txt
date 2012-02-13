# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (ApplicationIT_regularscalar) 

# set ITL path

# find the mpi package
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
# MOD-BY-LEETEN 02/13/2012-FROM:	# include_directories(MPI_INCLUDE_PATH)
include_directories(${MPI_INCLUDE_PATH})
if( DEFINED WIN32 )
	link_directories(${MPIW32_DIR}/lib)
endif()
# MOD-BY-LEETEN 02/13/2012-END
# use mpicxx to replace 
# DEL-BY-LEETEN 02/13/2012-BEGIN
	# SET(CMAKE_C_COMPILER mpicc)
	# SET(CMAKE_CXX_COMPILER mpicxx)
# DEL-BY-LEETEN 02/13/2012-END

# specify the path to the headers
include_directories(
		. 
		./include
		${ITL_DIR}/include
		# ${ITL_DIR}/include/datastore
		# ${ITL_DIR}/include/field
		# ${ITL_DIR}/include/grid
		# ${ITL_DIR}/include/metrics
		# ${ITL_DIR}/include/util
	)

link_directories(
	${ITL_DIR}/lib
)

# grab all the source codes
AUX_SOURCE_DIRECTORY ("src" 	SRC_FILES)

# build the project as a library
add_executable(${PROJECT_NAME} ${SRC_FILES})

# ADD-BY-LEETEN 02/13/2012-BEGIN
if( DEFINED WIN32 )	
target_link_libraries(${PROJECT_NAME} 
	mpi cxx
	)
endif()	
# ADD-BY-LEETEN 02/13/2012-END

# link the libraries
target_link_libraries(${PROJECT_NAME} ITLib) 

#
# $Log$
#

