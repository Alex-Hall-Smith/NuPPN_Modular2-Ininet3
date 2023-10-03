# Compiler and options for Intel
FC=ifort 

# ^_^ tell the code we're using ifort
CPPFLAGS=-DIFORT

FFLAGSOPT= -O3 
FFLAGSDEBUG=-g  -check all -warn all

ifdef LAPACK_LIBS
CPPFLAGS += -DLAPACK_LIBS
LDLIBS += $(LAPACK_LIBS)
ifdef LD_PATH_LAPACK
LDFLAGS += -L$(LD_PATH_LAPACK) -Wl,-rpath=$(LD_PATH_LAPACK)
endif
endif

