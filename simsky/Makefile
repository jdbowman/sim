# Makefile for simsky. Created Apr 2009 by JDB
# 2013-07-26:
# Add variable for HEALPix 2.15a dependent
# You must build the c++ library before compiling simsky
HPX_CXX_INC = ./Healpix_2.15a/src/cxx/generic_gcc/include
HPX_CXX_LIB = ./Healpix_2.15a/src/cxx/generic_gcc/lib

CXX = g++
CXX_INCS=-I$(HPX_CXX_INC)
CXXFLAGS = -g -Wall $(CXX_INCS) -L$(HPX_CXX_LIB)

all: srcgen galfix galmask	

srcgen: srcgen.cc
	$(CXX) $(CXXFLAGS) -o srcgen srcgen.cc -lgsl -lgslcblas -lm -lhealpix_cxx -lcxxsupport -lfftpack -lcfitsio

galfix: galfix.cc
	$(CXX) $(CXXFLAGS) -o galfix galfix.cc -lgsl -lgslcblas -lm -lhealpix_cxx -lcxxsupport -lfftpack -lcfitsio

galmask: galmask.cc
	$(CXX) $(CXXFLAGS) -o galmask galmask.cc -lgsl -lgslcblas -lm -lhealpix_cxx -lcxxsupport -lfftpack -lcfitsio

clean:
	rm -f *.o srcgen galfix galmask
