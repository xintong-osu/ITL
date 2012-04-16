CC = mpicxx
CCFLAGS = -O0 -Wall -g
INCLUDES += -I/home/abon/build/diy/include
#LIBDIRS += 
#LIBS += 

SOURCES = ./Main_reader.cpp
#      Main_reader.cpp

OBJECTS = Main_reader.o   	

# --- Make targets ---

default: all

all: serialio_diy

serialio_diy: $(OBJECTS)
	$(CC) $(LIBDIRS) $(INCLUDES) $(CCFLAGS) -o $@ $< $(LIBS)
	rm -f $(OBJECTS)

%.o: ./src/%.cpp
	$(CC) $(INCLUDES) $(CCFLAGS) -c $< -o $@ 

# --- Remove binary and executable files ---

clean:
	rm -f serialio_diy $(OBJECTS)

