AYB README File
===============

Prerequisites
-------------
The following utilities and libraries must be installed to make/run the program:

Utilities::
- gcc (GNU C compiler; part of GNU Compiler Collection)
- gfortran (GNU Fortran compiler; part of GNU Compiler Collection)
  (not required from v2.05)
- make (GNU make utility)

Libraries::
- BLAS (Basic Linear Algebra Subprograms)
- LAPACK (Linear Algebra PACKage)
- zlib (general purpose compression library)
- bzip2 (alternative compression library)

External Links::
- http://gcc.gnu.org/[GCC]
- http://www.gnu.org/software/make/[GNU Make]
- http://www.netlib.org/blas/[BLAS]
- http://www.netlib.org/lapack/[LAPACK]
- http://zlib.net/[zlib]
- http://bzip.org/[bzip2]


Obtaining AYB
-------------
Source code is freely available to download from <http://www.ebi.ac.uk/goldman-srv/AYB/>


Build AYB
---------
unzip the archive:

---------------------------------------------
$ tar -xzf AYB_current.tgz
---------------------------------------------

Go to the src directory. The makefile is suitable for both generic Linux and Mac Os-X.

---------------------------------------------
$ make
---------------------------------------------

The AYB executable will be in the bin directory.

{nbsp}

Notes::
If OpenMP libraries are not available then the dependency can be removed by 
deleting the -fopenmp option from the makefile.
