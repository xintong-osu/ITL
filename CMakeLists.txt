# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITL)

# ADD-BY-ABON 11/07/2011-BEGIN
set(WITH_ITL_API OFF CACHE BOOL "whether ITL API is compiled") 
set(BUILD_DIY_EXAMPLES OFF CACHE BOOL "whether DIY examples are to be built") 
# ADD-BY-ABON 11/07/2011-END
# ADD-BY-ABON 07/01/2012-BEGIN
set(BUILD_MOAB_EXAMPLES OFF CACHE BOOL "whether MOAB examples are to be built") 
# set(BUILD_COMPRESSION_EXAMPLES OFF CACHE BOOL "whether Compression examples are to be built") 
# ADD-BY-ABON 07/01/2012-END

# ADD-BY-LEETEN 08/15/2011-BEGIN
set(MYLIB_DIR 	"/homes/leeten/mylib" CACHE PATH "The path to MYLIB")
set(WITH_PNETCDF OFF CACHE BOOL "whether PNETCDF is used")
set(PNETCDF_DIR  "/usr/local/parallel-netcdf" CACHE PATH "installed path of PNETCDF")

# ADD-BY-LEETEN 09/01/2011-BEGIN
set(PNETCDF_LIB	${PNETCDF_DIR}/lib)
set(PNETCDF_INC	${PNETCDF_DIR}/include)

set(NETCDF_DIR  "/usr/local/netcdf" CACHE PATH "installed path of NETCDF")
set(NETCDF_INC  "${NETCDF_DIR}/include")
set(NETCDF_LIB  "${NETCDF_DIR}/lib")
# ADD-BY-LEETEN 09/01/2011-END

# ADD-BY-LEETEN 08/29/2011-BEGIN
set(HDF5_DIR 	"/usr/local/hdf5" CACHE PATH "The path to HDF5")
# ADD-BY-ABON 11/07/2011-BEGIN
set(DIY_DIR 	"/usr/local/diy" CACHE PATH "The path to DIY")

set(WITH_ZLIB   OFF CACHE BOOL "whether ZLIB is used") # ADD-BY-ABON 11/07/2011-BEGIN
set(ZLIB_DIR 	"/usr/local/zlib" CACHE PATH "The path to zlib")
# ADD-BY-ABON 11/07/2011-END

# Set Matlab engine
set( MATENGINE_INC_DIR "/home/abon/install/matlabr2011b/extern/include" CACHE PATH "Path to Matlab engine" )
set( MATENGINE_LIB_DIR "/home/abon/install/matlabr2011b/bin/glnxa64" CACHE PATH "Path to Matlab engine libs" )

# set the path of DIY files
set(DIY_INC "${DIY_DIR}/include")  
set(DIY_LIB "${DIY_DIR}/lib")  

# ADD-BY-ABON 11/07/2011-BEGIN
# set the path of BIL files
set(BIL_INC "${DIY_DIR}/src/bil-0.6.0/src" CACHE PATH "The path to bil headers" )  
set(BIL_LIB "${DIY_DIR}/src/bil-0.6.0" CACHE PATH "The path to bil library" )
# ADD-BY-ABON 11/07/2011-END

# ADD-BY-ABON 11/07/2011-BEGIN
# set the path to zlib
set(ZLIB_INC "${ZLIB_DIR}/include")  
set(ZLIB_LIB "${ZLIB_DIR}/lib" )
# ADD-BY-ABON 11/07/2011-END

# ADD-BY-ABON 07/01/2011-BEGIN
if( ${BUILD_MOAB_EXAMPLES} )

  set(MOAB_DIR "/homes/chaudhua/install" CACHE PATH "The installed path to MOAB" )
  set(MOAB_INC "${MOAB_DIR}/include")
  set(MOAB_LIB "${MOAB_DIR}/lib")

  set(VTK_DIR "/homes/chaudhua/install" CACHE PATH "The installed path to VTK" )
  set(VTK_INC "${VTK_DIR}/include/vtk-5.10")
  set(VTK_LIB "${VTK_DIR}/lib/vtk-5.10")

endif()

#if( ${BUILD_COMPRESSION_EXAMPLES} )
#  set( GSL_INC "/usr/include" CACHE PATH "The path to gsl headers" )  
#  set( GSL_LIB "usr/lib" CACHE PATH "The path to gsl library" )	
  
#  set( OPENCV_DIR "/homes/chaudhua/install" CACHE PATH "The installed path of OpenCV" )
#  set( OPENCV_INC "${OPENCV_DIR}/include" )
#  set( OPENCV_INC2 "${OPENCV_DIR}/include/opencv" )
#  set( OPENCV_INC3 "${OPENCV_DIR}/include/opencv2" )
#  set( OPENCV_LIB "${OPENCV_DIR}/lib")
  
#  set( MY_INSTALL_DIR "/homes/chaudhua/install" CACHE PATH "The installed path of my library" )
  
#endif()
# ADD-BY-ABON 07/01/2011-END

# set the path of ITL files
set(ITL_DIR "${CMAKE_SOURCE_DIR}" )
set(ITL_INC "${ITL_DIR}/include")  
set(ITL_LIB "${ITL_DIR}/lib")  
# ADD-BY-LEETEN 08/29/2011-END

# ADD-BY-LEETEN 02/13/2012-BEGIN
if( DEFINED WIN32 )
	set(MPIW32_DIR 	"C:/Program Files (x86)/DeinoMPI" CACHE PATH "The path to MPI on Windows 32 platform.")
endif()

add_definitions(
	-DMPICH_IGNORE_CXX_SEEK 
	-D_USE_MATH_DEFINES
	)
# ADD-BY-LEETEN 02/13/2012-END

# if(WITH_PNETCDF EQUAL ON)
if( ${WITH_PNETCDF} )
  add_definitions(-DWITH_PNETCDF)
  # specify the path to the headers
    include_directories(
      ${PNETCDF_INC}
    )
    set(NETCDF_LFLAGS "-L${PNETCDF_LIB} -lpnetcdf")
  else()
    # specify the path to the headers
    include_directories(
      ${NETCDF_INC}
      )
    set(NETCDF_LFLAGS "-L${NETCDF_LIB} -lnetcdf")
endif()
# ADD-BY-LEETEN 08/15/2011-END

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

add_subdirectory (include) # ADD-BY-LEETEN 09/07/2011 

add_subdirectory (examples/ApplicationIT_regularscalar)
add_subdirectory (examples/ApplicationIT_regularvector)
add_subdirectory (examples/ApplicationIT_irregularScalar)
# ADD-BY-ABON 11/07/2011-BEGIN
if( ${BUILD_DIY_EXAMPLES} )
	add_subdirectory (examples/ApplicationIT_diy_global)	
	add_subdirectory (examples/ApplicationIT_diy_blockwise)
	add_subdirectory (examples/ApplicationIT_diy_local)
	# DEL-BY-LEEETEN 01/04/2012: add_subdirectory (examples/ApplicationIT_diy_local_communication)
	add_subdirectory (examples/ApplicationIT_diy_crossvalidate)	
endif()
# ADD-BY-ABON 11/07/2011-END

# ADD-BY-LEETEN 08/29/2011-BEGIN
if( ${WITH_ITL_API} ) # MOD-BY-LEETEN 11/10/2011-FROM:	if( $WITH_ITL_API ) 
  add_subdirectory (examples/nek5000)
  add_subdirectory (examples/flash)
endif() 			
# ADD-BY-LEETEN 08/29/2011-END

# ADD-BY-ABON 07/01/2012-BEGIN
if( ${BUILD_MOAB_EXAMPLES} )
  add_subdirectory ( examples/ApplicationIT_nekmoab )	
endif()

#if( ${BUILD_COMPRESSION_EXAMPLES} )
#  add_subdirectory ( examples/ApplicationIT_compressionindexing )	
#endif()
# ADD-BY-ABON 07/01/2012-END

#
# $Log$
#

