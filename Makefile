QUIJOTEFCOMP=ifort
QUIJOTEINC=-I/nas/scratch/sharper/QUIJOTE/Pipeline/Mapping/lib/include -I/nas/scratch/sharper/QUIJOTE/Pipeline/Mapping/lib/include
QUIJOTELIB=-L/nas/scratch/sharper/QUIJOTE/Pipeline/Mapping/lib/lib -L/nas/scratch/sharper/QUIJOTE/Pipeline/Mapping/lib/lib -L/nas/scratch/sharper/QUIJOTE/Pipeline/DataAccess/lib

LIBS=-lsla -lhealpix -lcfitsio
FLAGS=

all: Ephem Coordinates Binning

Ephem:

	f2py -c --f90exec=$(QUIJOTEFCOMP) $(FLAGS) -m Ephem Ephem.f90 $(QUIJOTEINC) $(QUIJOTELIB) $(LIBS)

Coordinates:

	f2py -c --f90exec=$(QUIJOTEFCOMP) $(FLAGS) -m Coordinates_tpoint Coordinates_tpoint.f90 $(QUIJOTEINC) $(QUIJOTELIB) $(LIBS)

	f2py -c --f90exec=$(QUIJOTEFCOMP) $(FLAGS) -m Coordinates Coordinates.f90 $(QUIJOTEINC) $(QUIJOTELIB) $(LIBS)

Binning:

	f2py -c --f90exec=$(QUIJOTEFCOMP) $(FLAGS) -m Binning Binning.f90 $(QUIJOTEINC) $(QUIJOTELIB) $(LIBS)
