COMP=#Fortran Compiler command
INC=-I#/Path/To/Python/Include

all: Binning

Binning:

	f2py -c --f90exec=$(COMP) -m fBinning Tools/fBinning.f90  
	mv fBinning.so Tools/fBinning.so