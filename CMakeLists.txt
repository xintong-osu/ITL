# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (ITL)

# ADD-BY-LEETEN 08/15/2011-BEGIN
set(WITH_PNETCDF OFF CACHE BOOL "whether PNETCDF is used")
set(PNETCDF_PATH  "/usr/local/parallel-pnetcdf" CACHE PATH "installed path of PNETCDF")

if(WITH_PNETCDF EQUAL ON)
  add_definitions(-DWITH_PNETCDF)
  # specify the path to the headers
  include_directories(
    ${PNETCDF_PATH}/include
    )
endif()
# ADD-BY-LEETEN 08/15/2011-END

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

add_subdirectory (examples/ApplicationIT_regularvector)
add_subdirectory (examples/ApplicationIT_regularscalar)

#
# $Log$
#

