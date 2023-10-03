# Compiler and options for gfortran on Darwin
FC := gfortran

# ^_^ tell the code we're using gfort
CPPFLAGS=-DGFORT

WARN:= -Wall -Wtabs -Wno-unused-dummy-argument -Wno-unused-variable 
OPT:= -O3 -fforce-addr
DBG:= -g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
FFLAGSOPT := $(OPT)
FFLAGSDEBUG := $(WARN) $(DBG)


# ACML compatibility
ifdef LAPACK_LIBS
CPPFLAGS += -DLAPACK_LIBS
LDLIBS += $(LAPACK_LIBS)
ifdef LD_PATH_LAPACK
LDFLAGS += -L$(LD_PATH_LAPACK) -Wl,-rpath=$(LD_PATH_LAPACK)
endif
endif

