############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _linux

#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= $(SOLVER)/source/SuperLU_5.0
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_5.0.a

TMGLIB       	= libtmglib.a
LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc
CFLAGS       = -O3 -g
NOOPTS       = 
FORTRAN	     = gfortran #g77
FFLAGS       = -O2 -g # -fopenmp
LOADER       = $(CC)
LOADOPTS     = # -fopenmp

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/matlab

