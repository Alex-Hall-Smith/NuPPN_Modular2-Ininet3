# Compiler and options for gfortran on Linux 64 bit
FC := mpif90

# ^_^ tell the code we're using gfort
CPPFLAGS=-DGFORT

WARN:= -Wall -Wtabs -Wno-unused-dummy-argument -Wno-unused-variable 
OPT:= -O3 -fforce-addr
DBG:= -O0 -g -fbacktrace -fcheck=all $(WARN)
FFLAGSOPT := $(OPT)
FFLAGSDEBUG := $(DBG)


# ACML compatibility
ifdef LAPACK_LIBS
CPPFLAGS += -DLAPACK_LIBS
LDLIBS += $(LAPACK_LIBS)
ifdef LD_PATH_LAPACK
LDFLAGS += -L$(LD_PATH_LAPACK) -Wl,-rpath=$(LD_PATH_LAPACK)
endif
endif

LDFLAGS += -Wl,-rpath=$(SEHOME)/lib
