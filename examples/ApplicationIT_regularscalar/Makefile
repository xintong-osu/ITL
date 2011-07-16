include ../../user_defs.mk

# --- Definition of macros ---

# Franklin
#INCLUDES += -I../../headers -I../../headers/datastore -I../../headers/field \
#	    -I../../headers/grid -I../../headers/metrics -I../../headers/util \
#           -I./headers -I../ITLib/headers \
#LIBDIRS += -L../..
#LIBS += -lITLib 

# Davinci
CC = mpicxx
CCFLAGS = -O0 -g
INCLUDES += -I../../include -I./include
LIBDIRS += -L../../lib
LIBS += -lITLib -lmpich -lmpichcxx

# Exavis
#CC = mpicxx
#CCFLAGS = -O0 -g
#INCLUDES += -I../../headers -I../../headers/datastore -I../../headers/field \
	    -I../../headers/grid -I../../headers/metrics -I../../headers/util \
	    -I./headers
#LIBDIRS += -L../..
#LIBS += -lITLib -lmpich -lmpichcxx

SOURCES = ./src/*.cpp
#      MainIT_regscalar.cpp

OBJECTS = MainIT_regscalar.o   	

# --- Make targets ---

default: all

all: ApplicationIT_regularscalar

ApplicationIT_regularscalar: $(OBJECTS)
	$(CC) $(LIBDIRS) $(INCLUDES) $(CCFLAGS) -o $@ $< $(LIBS)
	rm -f $(OBJECTS)

%.o: ./src/%.cpp
	$(CC) $(INCLUDES) $(CCFLAGS) -c $< -o $@ 

# --- Remove binary and executable files ---

clean:
	rm -f AplicationIT_regularscalar $(OBJECTS)
	
	
