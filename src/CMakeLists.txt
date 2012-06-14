# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITLib) 

# set the path of MYLIB files
set(MYLIB_INC '${MYLIB_DIR}/include')  
set(MYLIB_LIB '${MYLIB_DIR}/lib')  

# find the mpi package
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
# MOD-BY-LEETEN 02/13/2012-FROM:	include_directories(MPI_INCLUDE_PATH)
include_directories(${MPI_INCLUDE_PATH})
# MOD-BY-LEETEN 02/13/2012-END

# use mpicxx to replace 
# DEL-BY-LEETEN 02/13/2012-BEGIN
	# SET(CMAKE_C_COMPILER mpicc)
	# SET(CMAKE_CXX_COMPILER mpicxx)
# DEL-BY-LEETEN 02/13/2012-END

# specify the path to the headers
include_directories(
		. 
		${ITL_INC}
		${MYLIB_INC}
	)

# grab all the source codes
list(APPEND SRC_FILES
	ITL_trianglepatch.cpp
	ITL_entropycore.cpp   
	ITL_vectormatrix.cpp
	ITL_histogram.cpp
	ITL_base.cpp         
	ITL_spacetreenode.cpp
	ITL_SphereSpace.cpp
)

# ADD-BY-LEETEN 02/13/2012-BEGIN
file(GLOB HEADER_FILES "../include/*.h" )
list(APPEND SRC_FILES ${HEADER_FILES})
# ADD-BY-LEETEN 02/13/2012-END

if( ${WITH_ITL_API} )
	list(APPEND SRC_FILES
		itl.cpp              
		ITL_random_field.cpp
	)
endif()



# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES})

# ADD-BY-LEETEN 02/13/2012-BEGIN
set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}_r"
)
# ADD-BY-LEETEN 02/13/2012-END

install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)
	
#
# $Log$
#

