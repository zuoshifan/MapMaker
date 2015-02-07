COMP=ifort
INC=-I/local/sharper/anaconda/include
FLAGS=

all: Binning

Binning:

	f2py -c --f90exec=$(COMP) $(FLAGS) -m fBinning fBinning.f90  
