# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (ApplicationIT_compressionindexing) 

# set ITL path
# DEL-BY-LEETEN 08/30/2011		set(ITL_DIR "../.." CACHE PATH "The path to ITL")

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
# ADD-BY-ABON 07/01/2012-BEGIN
include_directories( 
		${GSL_INC}
	)
# ADD-BY-ABON 07/01/2012-END

link_directories(
	${ITL_DIR}/lib
)
link_directories(
	${BIL_LIB}
)
if( ${WITH_ZLIB} ) # ADD-BY-LEETEN 01/04/2012
link_directories(
	${ZLIB_LIB}
)
# ADD-BY-LEETEN 01/04/2012-BEGIN
endif()

if( ${WITH_PNETCDF} )
# ADD-BY-LEETEN 01/04/2012-END
link_directories(
	${PNETCDF_LIB}
)
endif() # ADD-BY-LEETEN 01/04/2012D
link_directories(
	${DIY_LIB}
)
# ADD-BY-ABON 07/01/2012-BEGIN
link_directories(
	${GSL_LIB}
)
# ADD-BY-ABON 07/01/2012-END


# grab all the source codes
AUX_SOURCE_DIRECTORY ("src" 	SRC_FILES)

# build the project as a library
add_executable(${PROJECT_NAME} ${SRC_FILES})

# link the libraries
target_link_libraries(${PROJECT_NAME} ITLib) 
target_link_libraries(${PROJECT_NAME} diy) 
target_link_libraries(${PROJECT_NAME} bil)
target_link_libraries(${PROJECT_NAME} z) 
target_link_libraries(${PROJECT_NAME} pnetcdf) 
#
# $Log$
#

