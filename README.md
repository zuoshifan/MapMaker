MapMaker
========

Destriping Map-Making Software for astronomical data processing.

=======

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
