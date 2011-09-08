# Exavis
CC = g++
CCFLAGS = -O0 -Wall -g
INCLUDES = -I. -I/home/chaudhua/install/diy/diy/include -I/home/chaudhua/build/diy/include  
LIBDIRS = -L/homes/chaudhua/Install/diy/diy/lib 
LIBS = -ldiy

SOURCES = ./Main_reader.cpp
# Main_reader.cpp

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
