include ../../user_defs.mk

# --- Definition of macros ---

# Franklin
#INCLUDES += -I../../headers -I../../headers/datastore -I../../headers/field \
#	    -I../../headers/grid -I../../headers/metrics -I../../headers/util \
#           -I./headers -I../ITLib/headers \
#LIBDIRS += -L../..
#LIBS += -lITLib 

# Davinci
#CC = mpicxx
#CCFLAGS = -O0 -g
#INCLUDES += -I../../include -I./include
#LIBDIRS += -L../../lib
#LIBS += -lITLib -lmpich -lmpichcxx

# Exavis
CC = mpicxx
CCFLAGS = -O0 -g
INCLUDES += -I../../include -I./include -I/home/chaudhua/install/diy/diy/include -I/home/chaudhua/build/diy/include -I/home/chaudhua/build/diy/src/bil-0.6.0/src -I/home/chaudhua/install/zlib-1.2.5/include
LIBDIRS += -L../../lib -L/home/chaudhua/install/diy/diy/lib -L/home/chaudhua/build/diy/src/bil-0.6.0 -L/home/chaudhua/install/zlib-1.2.5/lib
LIBS += -lITLib -ldiy -lbil -lz

# Abon's desktop at ANL
#INCLUDES += -I../../include -I./include -I/homes/chaudhua/Install/diy/diy/include -I/homes/chaudhua/Build/diy_build/include -I/homes/chaudhua/Build/diy_build/src/bil-0.6.0/src -I/homes/chaudhua/Install/zlib-1.2.5/include
#LIBDIRS += -L../../lib -L/homes/chaudhua/Install/diy/diy/lib -L/homes/chaudhua/Build/diy_build/src/bil-0.6.0 -L/homes/chaudhua/Install/zlib-1.2.5/lib
#LIBS += -lITLib -ldiy -lbil -lz
#LIBS += -lITLib -ldiy -lbil

SOURCES = ./src/MainIT_diy_local_communication.cpp
#      MainIT_diy_local_communication.cpp

OBJECTS = MainIT_diy_local_communication.o   	

# --- Make targets ---

default: all

all: ApplicationIT_diy_local_communication

ApplicationIT_diy_local_communication: $(OBJECTS)
	$(CC) $(LIBDIRS) $(INCLUDES) $(CCFLAGS) -o $@ $< $(LIBS)
	rm -f $(OBJECTS)

%.o: ./src/%.cpp
	$(CC) $(INCLUDES) $(CCFLAGS) -c $< -o $@ 

# --- Remove binary and executable files ---

clean:
	rm -f ApplicationIT_diy_local_communication $(OBJECTS)
	
	
