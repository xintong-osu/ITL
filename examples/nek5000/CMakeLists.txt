# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (itl-nek5000-examples) 

# DEL-BY-LEETEN 09/01/2011-BEGIN
# if ( ${WITH_PNETCDF} )
#	SET(NETCDF_LFLAGS	"-L${PNETCDF_LIB} -lpnetcdf")
# else()
#	SET(NETCDF_LFLAGS	"-lnetcdf")
# endif()
# DEL-BY-LEETEN 09/01/2011-END

SET(ITL_LFLAGS "\"-L${ITL_LIB} -lITLib -L${DIY_LIB} -ldiy -lstdc++ ${NETCDF_LFLAGS}\"")

# configure the NEK5000 makefile 'makenek'
configure_file (
  "makenek.in"
  "makenek"
  )
 