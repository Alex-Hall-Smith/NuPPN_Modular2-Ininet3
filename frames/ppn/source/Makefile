# Copy Make.local.sample to Make.local, for example if you are on
# helix there may already be a Make.local.helix file, in which case
# you type
# $ cp Make.local.helix Make.local
# and edit to suit your local configuration.
# You don't have to 'make' the folder anymore.

MAIN := ppn
OBJS := frame_knobs.o output.o ppn.o

-include Make.local
PHYSICS = $(PPN)/physics
SOLVER = $(PPN)/solver
UTILS = $(PPN)/utils

ifndef ARCH
$(error ARCH is not defined; should be defined in Make.local)
endif
include ./ARCH/$(ARCH)/Makefile

CWD := $(shell pwd)
ONEUPWD := $(shell (cd ..;pwd))

SUBDIRS = $(UTILS)/source $(PHYSICS)/source $(SOLVER)/source

SOLVERDEPS:= frame_knobs.o frame_knobs.mod
PHYSICSDEPS:= frame_knobs.o frame_knobs.mod

LDLIBS += -lsolver -lphysics -lutils
FFLAGS += $(foreach dir,$(SUBDIRS),-I$(dir)) -I$(CWD) 
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

# let the code know this is PPN
CPPFLAGS+=-DPPN

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

ppn.o: ppn.F90
	$(FC) $(FFLAGS) -c $^

output.mod: output.F90 output.o
	$(FC) $(FFLAGS) -c $<

output.o: output.F90
	$(FC) $(FFLAGS) -c $^

frame_knobs.mod: frame_knobs.F90 frame_knobs.o
	$(FC) $(FFLAGS) -c $<

frame_knobs.o: frame_knobs.F90
	$(FC) $(FFLAGS) -c $^

$(UTILS)/source:
	$(MAKE) -C $@ PPN='true'

$(PHYSICS)/source: $(PHYSICSDEPS)
	$(MAKE) -C $@
	rsync -v --ignore-existing $(PHYSICS)/run_template/* ../run_template
	rsync -v --ignore-existing $(SOLVER)/run_template/* ../run_template
	rm -f  $(ONEUPWD)/NPDATA
	ln -s  $(PHYSICS)/NPDATA  $(ONEUPWD)/NPDATA

$(SOLVER)/source: $(SOLVERDEPS)
	$(MAKE) -C $@
	rsync -v --ignore-existing $(SOLVER)/run_template/* ../run_template

$(SUPERLU_PATH)/FORTRAN/c_fortran_dgssv.o:
	$(info now building superLU FORTRAN interface ... )
	$(MAKE) -C $(SUPERLU_PATH)/FORTRAN/ c_fortran_dgssv.o

$(SUPERLU_PATH)/lib/libsuperlu_5.0.a:
	$(info superLU libraries did not exist... trying to build them ... )
	$(MAKE) BLASLIB=$(BLASLIB) SOLVER=$(SOLVER) -C $(SUPERLU_PATH) superlulib

$(SOLVER)/source/slu_f90/slu.o: $(SOLVER)/source/slu_f90/slu.f90
	$(FC) $(FFLAGS) -c -o $@ $^

parameter.inc :
	cat $(PHYSICS)/source/parameter_physics.inc ./parameter_frame.inc > ./parameter.inc

clean:
	-$(RM) *.o $(MAIN).exe *pyppn* *% *~ *genmod.mod *genmod.f90

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
