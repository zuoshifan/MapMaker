MapMaker
========

Map-making package for single dish astronomical data.


Installation
============

Prerequistes
------------

The code in this repository was written for Python 2.7 and depends on
the following Python modules:

    NumPy
    SciPy

    mpi4py (optional)

Installation
------------

Change to the top directory of the package, install it by the usual
methods, either the standard ::

    $ python setup.py install [--user]

or in the develop mode ::

    $ python setup.py develop [--user]

It should also be installable directly with `pip` using the command ::

    $ pip install [-e] git+https://github.com/zuoshifan/MapMaker.git

Testing installation
--------------------

Open Python or iPython in a terminal and type

    $ import MapMaker

If everything has installed correctly this should work with no errors.


Change-Log
==========

V0.4: Changes:

    - Re-organize it to be a installable package.

V0.3: Changes:

	- Re-organised structure so that Destriper and ML mappers both share
      the same tools and CGM codes.

	- Added tutorial and instructions to wiki.

V0.2: Changes:

	- Renamed Binning.py to nBinning.py and CGM.py to nCGM.py.

	- Fixed a bug where the hitmap was not being calculated.

	- Destriping is fully functional on a single thread. MPI
	  functions not yet tested.

V0.1a: Changes:

	- Fixed a number of problems with the code. Everything should
	  at least run correctly now. No guarantees it will actually
	  destripe anything though.

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

V0.0: Initial upload of the code. Not exactly sure if all these
      versions work well together. E.g. This map-maker probably does not
      work as it is.


