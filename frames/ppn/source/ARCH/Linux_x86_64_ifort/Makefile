# Compiler and options for Intel
FC = $(MPIHOME)/bin/mpif90

# ^_^ tell code that we're using ifort
CPPFLAGS=-DIFORT

FFLAGSOPT = -O3  
FFLAGSDEBUG =-g  -check all -warn all

ifdef LAPACK_LIBS
CPPFLAGS += -DLAPACK_LIBS
LDLIBS += $(LAPACK_LIBS)
ifdef LD_PATH_LAPACK
LDFLAGS += -L$(LD_PATH_LAPACK) -Wl,-rpath=$(LD_PATH_LAPACK)
endif
endif

# Support large networks: allow data to be stored outside the first
# 2GB of address space.
FFLAGS += -mcmodel=medium -shared-intel
LDFLAGS += -mcmodel=medium -shared-intel
