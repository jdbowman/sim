Simsky generate foreground model suitable for further low-frequency radio
astronomical modeling. It consists of 3 programs: 1) srcgen: generate point
sources as well as source confusion map. 2) galmask and 3) galfix can be used
together to produce galactic diffuse emission.

Note: According to Judd, galmask and galfix do not work very well, and they
have not been tested by me (Piyanat).

Original author: Judd Bowman
Last edited: 2013-07-18 by Piyanat Kittiwisit

INSTALLATION
------------
All program required HEALPix and CFITSIO packages. The most recent version of
the two packages that work with simsky (HEALPix 2.15a, cfitsio 3.2.40) are
included in the main directory.

To build the programs, you must first compile HEALPix c++ library and edit
the Makefile in the main directory to point to the HEALPix library.

KNOWN BUGS
----------
- Source clustering algorithm in srcgen is not working vert well. The program
produce an output map with recurring clustering pattern