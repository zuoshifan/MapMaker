COMP=gfortran
INC=-I/PATH/TO/python.h
FLAGS=

all: Binning

Binning:

	f2py -c --f90exec=$(FCOMP) $(FLAGS) -m fBinning fBinning.f90 $(INC) $(LIB)
