MapMaker
========

Destriping Map-Making Software for astronomical data processing.

=======

Installation (See wiki for example use)

Prerequistes

The code in this repository was written for Python 2.7.3 and depends on the following Python modules:

    NumPy
    SciPy

    mpi4py

    f2py (to compile Fortran modules)

If at all possible I highly recommend installing the Anaconda Python package, which will install all these modules and f2py by default.
Install Procedure

First the f2py Fortran shared libraries must be compiled.

    Enter MapMaker/
    Open: Makefile
    Edit: COMP= to equal the fortran compiler command for your system
    Edit: INC= to equal the directory of your Python include file
    Save and close Makefile
    Type: make

This should compile the Fortran libraries. Next you must add the directory above MapMaker to your PYTHONPATH

    If using bash type: export PYTHONPATH=$PYTHONPATH:/path/to/MapMaker
    If using csh type: setenv PYTHONPATH /path/to/MapMaker:$PYTHONPATH

Testing installation:

Open Python or iPython in a terminal and type: import MapMaker . If everything has installed correctly this should work with no errors.


======= Change-Log

V0.0: Initial upload of the code. Not exactly sure if all these
versions work well together. E.g. This map-maker probably does not
work as it is.

V0.1: All the correct code has now been uploaded. Changes:
 
      - Removed alot of excess code.
 
      - Added a class to hold the output maps and auxiliaries.

      - CGM function now expects two user written functions to derive
      Ax = b. One to define Ax and another to define b.  CGM also
      expects an arguments list that is passed to the two functions.

      - Created a new binning routine that interfaces with fBinning
      Fortran library. Most maps should be modified 'inplace'
      therefore minimising the memory footprint.

      - fBinning modified to contain only functions necessary for the
      map-maker. Several functions renamed. Made it possible to pass
      an array shorter than TOD to bin_map_ext to produce baseline
      maps.

      NB: This version probably still does not work. Everything is
      untested.

V0.1a: Changes:
	
	-Fixed a number of problems with the code. Everything should
	at least run correctly now. No guarantees it will actually
	destripe anything though.

V0.2: Changes:

	-Renamed Binning.py to nBinning.py and CGM.py to nCGM.py.

	-Fixed a bug where the hitmap was not being calculated.

	-Destriping is fully functional on a single thread. MPI
	functions not yet tested.
	
V0.3: Changes:
	
	-Re-organised structure so that Destriper and ML mappers both share the same tools and CGM codes.
	-Added tutorial and instructions to wiki.
