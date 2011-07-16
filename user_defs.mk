#----------------------------------------------------------------------------
#
# makefile user definitions
#
#----------------------------------------------------------------------------
#
# users: set paths in this file only
# you should not need to touch the other makefiles in the project
#
#----------------------------------------------------------------------------

# 1. Set your paths here 
# keep a copy of them below for later reuse, commenting others' out
# please don't delete other people's settings

# Franklin
#CC =CC
#CCFLAGS = -O0 -g
#LIBDIRS += 

# Exavis
#CC = mpicxx
#CCFLAGS = -O0 -g
#INCLUDES += -I/usr/include/mpich2  
#LIBDIRS += -L/usr/lib

# Abon's Linux
#CC = mpicxx
#CCFLAGS = -O0 -g
#INCLUDES += -I/home/abon/Install/mpich2-1.3.2p1/include
#LIBDIRS += -L/home/abon/Install/mpich2-1.3.2p1/lib

# Abon's Linux at ANL
CC = mpicxx
CCFLAGS = -O0 -Wall -g

