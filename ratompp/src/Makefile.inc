#
# Macros for g++ compiler 
#

# Path to LAPACK library
LAPACK    = $(HOME)/lib-gcc/lapack-3.2.2/liblapack.a
# Path to BLAS library
BLAS      = $(HOME)/lib-gcc/lapack-3.2.2/libblas.a

CXX       = g++
CXXFLAGS  = -DNDEBUG -O2 -funroll-all-loops -ffast-math -Wall 
# CXXFLAGS  =  -g -O -Wall -D_DEBUG
CXXINCL   = 
CXXLIB    = $(LAPACK) $(BLAS) -lgfortran

