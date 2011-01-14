#----------------------------------------------------------------------------
#
# Tom Peterka
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# Copyright Notice
# + 2010 University of Chicago
#
#----------------------------------------------------------------------------

INCLUDE = -I. -I../include
CCFLAGS = -g -c
C++ = mpicxx

OBJS = itl.o

.SUFFIXES: .cpp

default: all

all: libitl.a

Makefile: Makefile.in ../config.status
	../config.status $@

config.status: ../configure
	../config.status --recheck

$(OBJS): ../include/*.h

.cpp.o:
	$(C++) $(CCFLAGS) $(INCLUDE) $<

libitl.a : $(OBJS)
	rm -f $@
	ar cru $@ $(OBJS) 
	mv libitl.a ../lib
clean:
	rm -f *.o *.a
