include user_defs.mk

# --- Definition of macros ---
CCFLAGS = -O0 -g
INCLUDES += -I./include 

SOURCES = ./src/*.cpp
#      ITL_base.cpp ITL_entropycore.cpp\
#      ITL_vectormatrix.cpp ITL_histogram.cpp

OBJECTS = ITL_base.o ITL_vectormatrix.o ITL_histogram.o ITL_entropycore.o

# --- Make targets ---

default: all

all: ./lib/libITLib.a

./lib/libITLib.a: $(OBJECTS)
	ar cru $@ $(OBJECTS)
	rm -f $(OBJECTS)

%.o: ./src/%.cpp
	$(CC) $(INCLUDES) $(CCFLAGS) -c $< -o $@ 

# --- Remove binary and executable files ---

clean:
	rm -f ./lib/*.a $(OBJECTS)
	
	
