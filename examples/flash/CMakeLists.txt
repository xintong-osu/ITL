# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (itl-flash-examples) 

find_package(MPI REQUIRED)

if ( ${WITH_PNETCDF} )
	SET(NETCDF_LFLAGS	"-L${PNETCDF_LIB} -lpnetcdf")
else()
	SET(NETCDF_LFLAGS	"-lnetcdf")
endif()

SET(MPI_PATH ${MPI_DIR})
SET(HDF5_PATH ${HDF5_DIR})
SET(NCMPI_PATH ${PNETCDF_DIR})
SET(ITL_LFLAGS "-L${ITL_LIB} -lITLib -L${DIY_LIB} -ldiy -lstdc++ ${NETCDF_LFLAGS}")

# configure the makefile 
configure_file (
  "Makefile.h.in"
  "Makefile.h"
  )
 
