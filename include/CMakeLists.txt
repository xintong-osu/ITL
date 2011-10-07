cmake_minimum_required (VERSION 2.8)
project (ITLib-include) 

# grab all the source codes
file(GLOB HEADER_FILES "*.h" )

install(
	FILES ${HEADER_FILES}
	DESTINATION include)
	
#
# $Log$
#

