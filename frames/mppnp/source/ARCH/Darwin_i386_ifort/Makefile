# Compiler and options for Intel
FC = $(MPIHOME)/bin/mpif90

# ^_^ tell the code we're using ifort
CPPFLAGS=-DIFORT

FFLAGSOPT= -O3  -diag-disable vec
FFLAGSDEBUG=-g  -check all -warn all

ifdef LAPACK_LIBS
CPPFLAGS += -DLAPACK_LIBS
ifdef INCLUDE_LAPACK
FFLAGS += -I$(INCLUDE_LAPACK)
endif
LDLIBS += $(LAPACK_LIBS)
ifdef LD_PATH_LAPACK
LDFLAGS += -L$(LD_PATH_LAPACK)
endif
endif
