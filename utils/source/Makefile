OBJS = array_sizes.o \
       errors.o \
       sort.o \
       utils.o \
       constants.o \
       interp.o \
       akima.o \
       cubinterp.o \
       trajectories.o \
       communication.o

# hdf5 libraries
ifdef USE_HDF5
OBJS += fse_wrapper.o
endif

libutils.a: $(OBJS) $(patsubst %.o,libutils.a(%.o),$(OBJS))
	ranlib libutils.a

clean:
	$(RM) *.o *.mod libutils.a

%.o %.mod: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

