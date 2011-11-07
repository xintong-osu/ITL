# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITL)

set(WITH_ITL_API OFF CACHE BOOL "whether ITL API is compiled") # TEST-ADD
set(BUILD_DIY_EXAMPLES OFF CACHE BOOL "whether DIY examples are to be built") # TEST-ADD

# ADD-BY-LEETEN 08/15/2011-BEGIN
set(MYLIB_DIR 	"/homes/leeten/mylib" CACHE PATH "The path to MYLIB")
set(WITH_PNETCDF OFF CACHE BOOL "whether PNETCDF is used")
# MOD-BY-LEETEN 08/29/2011-FROM:	set(PNETCDF_DIR  "/usr/local/parallel-pnetcdf" CACHE PATH "installed path of PNETCDF")
# TO:
set(PNETCDF_DIR  "/usr/local/parallel-netcdf" CACHE PATH "installed path of PNETCDF")
# MOD-BY-LEETEN 08/29/2011-END	

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
set(ZLIB_DIR 	"/usr/local/zlib" CACHE PATH "The path to zlib")
# ADD-BY-ABON 11/07/2011-END

# set the path of DIY files
set(DIY_INC "${DIY_DIR}/include")  
set(DIY_LIB "${DIY_DIR}/lib")  

# ADD-BY-ABON 11/07/2011-BEGIN
# set the path of BIL files
set(BIL_INC "${DIY_DIR}/src/bil-0.6.0/src")  
set(BIL_LIB "${DIY_DIR}/src/bil-0.6.0" )
# ADD-BY-ABON 11/07/2011-END

# ADD-BY-ABON 11/07/2011-BEGIN
# set the path to zlib
set(ZLIB_INC "${ZLIB_DIR}/include")  
set(ZLIB_LIB "${ZLIB_DIR}/lib" )
# ADD-BY-ABON 11/07/2011-END

# set the path of ITL files
set(ITL_DIR "${CMAKE_SOURCE_DIR}")
set(ITL_INC "${ITL_DIR}/include")  
set(ITL_LIB "${ITL_DIR}/lib")  
# ADD-BY-LEETEN 08/29/2011-END

# if(WITH_PNETCDF EQUAL ON)
if( ${WITH_PNETCDF} )
  add_definitions(-DWITH_PNETCDF)
  # specify the path to the headers
  # MOD-BY-LEETEN 09/01/2011-FROM:
    # include_directories(
      # ${PNETCDF_DIR}/include
    # )
  # MOD-BY-LEETEN 09/01/2011-TO:
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
  # MOD-BY-LEETEN 09/01/2011-END
endif()
# ADD-BY-LEETEN 08/15/2011-END

# DEL-BY-LEETEN 09/01/2011-BEGIN
# # ADD-BY-LEETEN 08/29/2011-BEGIN
# set(PNETCDF_LIB	${PNETCDF_DIR}/lib)
# set(PNETCDF_INC	${PNETCDF_DIR}/include)
# # ADD-BY-LEETEN 08/29/2011-END
# DEL-BY-LEETEN 09/01/2011-END

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

add_subdirectory (include) # ADD-BY-LEETEN 09/07/2011 

add_subdirectory (examples/ApplicationIT_regularscalar)
add_subdirectory (examples/ApplicationIT_regularvector)
# ADD-BY-ABON 11/07/2011-BEGIN
if( ${BUILD_DIY_EXAMPLES} )
	add_subdirectory (examples/ApplicationIT_diy_global)	
	add_subdirectory (examples/ApplicationIT_diy_blockwise)
	add_subdirectory (examples/ApplicationIT_diy_local)
	add_subdirectory (examples/ApplicationIT_diy_local_communication)
endif()
# ADD-BY-ABON 11/07/2011-END

# ADD-BY-LEETEN 08/29/2011-BEGIN
if( $WITH_ITL_API ) # TEST-ADD
  add_subdirectory (examples/nek5000)
  add_subdirectory (examples/flash)
endif() 			# TEST-ADD
# ADD-BY-LEETEN 08/29/2011-END

#
# $Log$
#

