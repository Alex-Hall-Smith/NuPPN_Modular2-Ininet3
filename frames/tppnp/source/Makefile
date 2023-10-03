# Copy Make.local.sample to Make.local and edit to suit your local
# configuration, then run
#
# % make
#
# to build an optimized binary and libraries, or
#
# % make debug
#
# to build debug versions.

MAIN := tppnp
OBJS := traj_data.o frame_knobs.o tppnp.o

-include Make.local
PHYSICS = $(PPN)/physics
SOLVER  = $(PPN)/solver
UTILS   = $(PPN)/utils

ifndef ARCH
$(error ARCH is not defined; should be defined in Make.local)
endif
include ./ARCH/$(ARCH)/Makefile

CWD     := $(shell pwd)

SUBDIRS = $(UTILS)/source $(PHYSICS)/source $(SOLVER)/source

SOLVERDEPS  := frame_knobs.o frame_knobs.mod
PHYSICSDEPS := frame_knobs.o frame_knobs.mod

ifeq ($(USE_SE),YES)
LDLIBS   += -lse
FFLAGS   += -I$(SEHOME)/include
LDFLAGS  += -L$(SEHOME)/lib
CPPFLAGS += -DUSE_HDF5
endif

LDLIBS  += -lsolver -lphysics -lutils
FFLAGS  += $(foreach dir,$(SUBDIRS),-I$(dir)) -I$(CWD)
LDFLAGS += $(foreach dir,$(SUBDIRS),-L$(dir))

# superLU (open source sparse solver) compatibility
ifeq ($(USE_SUPERLU),YES)
SUPERLU_PATH = $(SOLVER)/source/SuperLU_5.0
CPPFLAGS+=-DUSE_SUPERLU
# for the superLU libraries and interfaces
OBJS += $(SUPERLU_PATH)/FORTRAN/c_fortran_dgssv.o
OBJS += $(SUPERLU_PATH)/lib/libsuperlu_5.0.a
LDLIBS += $(BLASLIB) -lm
# for the superLU interface by P. Edelmann
CFLAGS = -Wall -g -O3 -I$(SUPERLU_PATH)/SRC
OBJS += $(SOLVER)/source/slu_f90/slu_f90.o
OBJS += $(SOLVER)/source/slu_f90/slu.o
LDLIBS += -L$(SUPERLU_PATH)/lib -lsuperlu_5.0
SOLVERDEPS += $(SOLVER)/source/slu_f90/slu.o
SOLVERDEPS += $(SUPERLU_PATH)/FORTRAN/c_fortran_dgssv.o
endif

# tell the code that this is TPPNP
CPPFLAGS+=-DTPPNP

# pass NNN, NVCP and INDEX_REALIB to the preprocessor if the code
# is being compiled in the run directory
ifdef pNNN
ifdef pIDX_RCLB
ifdef pNVCP
CPPFLAGS+=-DpNNN=$(pNNN) -DpIDX_RCLB=$(pIDX_RCLB) -DpNVCP=$(pNVCP)
endif
endif
endif

export FC
export FFLAGS
export CPPFLAGS

.PHONY: all debug $(SUBDIRS)

all: FFLAGS += $(FFLAGSOPT)
all: $(MAIN).exe

debug: FFLAGS += $(FFLAGSDEBUG)
debug: $(MAIN).exe

$(MAIN).exe: parameter.inc $(SUBDIRS) $(OBJS) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

tppnp.o: tppnp.F90
	$(FC) $(FFLAGS) -c $^

frame_knobs.mod: frame_knobs.F90 frame_knobs.o
	$(FC) $(FFLAGS) -c $<

frame_knobs.o: frame_knobs.F90
	$(FC) $(FFLAGS) -c $^

traj_data.o: traj_data.F90
	$(FC) -cpp $(CPPFLAGS) $(FFLAGS) -c $^

$(UTILS)/source:
	$(MAKE) -C $@

$(PHYSICS)/source: $(PHYSICSDEPS)
	$(MAKE) -C $@

$(SOLVER)/source: $(SOLVERDEPS)
	$(MAKE) -C $@

$(SUPERLU_PATH)/FORTRAN/c_fortran_dgssv.o:
	$(info now building superLU FORTRAN interface ... )
	$(MAKE) -C $(SUPERLU_PATH)/FORTRAN/ c_fortran_dgssv.o

$(SUPERLU_PATH)/lib/libsuperlu_5.0.a:
	$(info superLU libraries did not exist... trying to build them ... )
	$(MAKE) BLASLIB=$(BLASLIB) SOLVER=$(SOLVER) -C $(SUPERLU_PATH) superlulib

$(SOLVER)/source/slu_f90/slu.o: $(SOLVER)/source/slu_f90/slu.f90
	$(FC) $(FFLAGS) -c -o $@ $^

parameter.inc:
	cat $(PHYSICS)/source/parameter_physics.inc ./parameter_frame.inc > ./parameter.inc

clean:
	-$(RM) *.o $(MAIN).exe *% *~ *genmod.mod *genmod.f90
	rm -f ./*.mod ./*.o

# superLU libraries clean

superclean:
	$(MAKE) SOLVER=$(SOLVER) BLASLIB=$(BLASLIB) -C $(SOLVER)/source/SuperLU_5.0 clean

distclean: clean
	$(MAKE) -C $(UTILS)/source clean
	$(MAKE) -C $(PHYSICS)/source clean
	$(MAKE) -C $(SOLVER)/source clean
	$(RM) $(SOLVER)/source/slu_f90/*.o
	$(RM) parameter.inc
	$(RM) ./slu.mod
