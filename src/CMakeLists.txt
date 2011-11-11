# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITLib) 

# MOD-BY-LEETEN 08/29/2011-FROM:
	# set(ITL_DIR "/homes/leeten/workspace/itl" CACHE PATH "The path to ITL")
	# set(MYLIB_DIR "/homes/leeten/mylib" CACHE PATH "The path to MYLIB")
# MOD-BY-LEETEN 08/29/2011-TO:
# set the path of MYLIB files
set(MYLIB_INC '${MYLIB_DIR}/include')  
set(MYLIB_LIB '${MYLIB_DIR}/lib')  
# MOD-BY-LEETEN 08/29/2011-END

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
		# MOD-BY-LEETEN 08/29/2011-FROM:	
			# ../include
			# ${ITL_DIR}/include
			# ${MYLIB_DIR}/include
		# MOD-BY-LEETEN 08/29/2011-TO:
		${ITL_INC}
		${MYLIB_INC}
		# MOD-BY-LEETEN 08/29/2011--END
	)

# grab all the source codes
# MOD-BY-LEETEN 11/10/2011-FROM:
	# AUX_SOURCE_DIRECTORY ("." 	SRC_FILES)
# TO:
list(APPEND SRC_FILES
	ITL_entropycore.cpp   
	ITL_vectormatrix.cpp
	ITL_histogram.cpp
	ITL_base.cpp         
	ITL_spacetreenode.cpp
	ITL_SphereSpace.cpp
)

if( ${WITH_ITL_API} )
	list(APPEND SRC_FILES
		itl.cpp              
		ITL_random_field.cpp
	)
endif()
# MOD-BY-LEETEN 11/10/2011-END



# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES})

# MOD-BY-LEETEN 09/07/2011-FROM:
	# # output the library to the lib/ folder
	# # MOD-BY-LEETEN 08/29/2011-FROM:	set(LIBRARY_OUTPUT_PATH "../lib")
	# # TO:
	# set(LIBRARY_OUTPUT_PATH ${ITL_LIB})
	# # MOD-BY-LEETEN 08/29/2011-END
# MOD-BY-LEETEN 09/07/2011-TO:
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)
	
# MOD-BY-LEETEN 09/07/2011-END

#
# $Log$
#

