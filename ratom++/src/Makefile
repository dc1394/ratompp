#
# Makefile for RAtom program
# Zbigniew Romanowski [romz@wp.pl]
#


include Makefile.inc

# All object files included into program
OBJ := $(patsubst %.cpp,%.o,$(wildcard ./Alg/*.cpp))    \
       $(patsubst %.cpp,%.o,$(wildcard ./ExCorr/*.cpp)) \
       $(patsubst %.cpp,%.o,$(wildcard ./Fem1D/*.cpp))  \
       $(patsubst %.cpp,%.o,$(wildcard ./Integ/*.cpp))  \
       $(patsubst %.cpp,%.o,$(wildcard ./Ks/*.cpp))     \
       $(patsubst %.cpp,%.o,$(wildcard ./Util/*.cpp))


# Path to final binary file (result of the compilation)
BINOUT = ../bin/RAtom.x


.PRECIOUS: $(BINOUT)
# .SILENT:

all:
	cd Alg;    make
	cd ExCorr; make
	cd Fem1D;  make
	cd Integ;  make
	cd Ks;     make
	cd Util;   make
	$(CXX) $(OBJ) $(CXXFLAGS) $(CXXLIB) -o $(BINOUT)


clean:
	rm -f $(OBJ) $(BINOUT)

