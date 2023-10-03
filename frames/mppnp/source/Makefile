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

MAIN := mppnp

OBJS := mixing.o fse_wrapper.o frame_knobs.o mppnp.o

SUBDIRS = $(UTILS)/source $(PHYSICS)/source $(SOLVER)/source
LDLIBS += -lsolver -lphysics -lutils -lse

-include Make.local
PHYSICS = $(PPN)/physics
SOLVER = $(PPN)/solver
UTILS = $(PPN)/utils
include ./ARCH/$(ARCH)/Makefile

SOLVERDEPS:= frame_knobs.o frame_knobs.mod
PHYSICSDEPS:= frame_knobs.o frame_knobs.mod

FFLAGS += $(foreach dir,$(SUBDIRS),-I$(dir)) -I$(CWD) -I$(SEHOME)/include -mcmodel=medium
LDFLAGS += $(foreach dir,$(SUBDIRS),-L$(dir)) -L$(SEHOME)/lib

CWD := $(shell pwd)
ONEUPWD := $(shell (cd ..;pwd))

ifndef ARCH
$(error ARCH is not defined; must be defined in Make.local)
endif
ifndef SEHOME
$(error SEHOME is not defined; must be defined in Make.local)
endif

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

# tell the code that this is MPPNP, so it knows to read
# the parameters regarding the multi-zone grid
CPPFLAGS+=-DMPPNP

# pass NNN, INDEX_REALIB, NVCP, msl, nrefmax and gfdim
#  to the preprocessor if the code
# is being compiled in the RUN directory
ifdef pNNN
ifdef pIDX_RCLB
ifdef pmsl
ifdef pgfdim
ifdef pnrefmax
ifdef pNVCP
CPPFLAGS+=-DpNNN=$(pNNN) -DpIDX_RCLB=$(pIDX_RCLB) -DpNVCP=$(pNVCP)
CPPFLAGS+=-Dpmsl=$(pmsl) -Dpgfdim=$(pgfdim) -Dpnrefmax=$(pnrefmax)
endif
endif
endif
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

mppnp.o: mppnp.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $^

mixing.mod: mixing.F90 mixing.o
	$(FC) $(FFLAGS) -c $<

mixing.o: mixing.F90
	$(FC) $(FFLAGS) -c $^

fse_wrapper.mod: fse_wraper.f90 fse_wrapper.o
	$(FC) --fixed-form -I$(SEHOME)/include -c $<

fse_wrapper.o: fse_wrapper.f90
	$(FC) --fixed-form -I$(SEHOME)/include -c $^

frame_knobs.mod: frame_knobs.F90 frame_knobs.o
	$(FC) $(FFLAGS) -c $<

frame_knobs.o: frame_knobs.F90
	$(FC) $(FFLAGS) -c $^

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
	-$(RM) *.o *.mod $(MAIN).exe *% *~ *genmod.mod *genmod.f90

# superLU libraries clean

superclean:
	$(MAKE) SOLVER=$(SOLVER) BLASLIB=$(BLASLIB) -C $(SOLVER)/source/SuperLU_5.0 clean

distclean: clean
	$(MAKE) -C $(UTILS)/source clean
	$(MAKE) -C $(PHYSICS)/source clean
	$(MAKE) -C $(SOLVER)/source clean
	rm -f $(SOLVER)/source/slu_f90/*.o
	rm -f parameter.inc
	rm -f ./slu.mod
