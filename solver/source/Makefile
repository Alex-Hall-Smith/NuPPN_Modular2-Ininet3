OBJS := lapack.o \
	solver_diagnostics.o \
	solver_knobs.o \
	sparsekit_module.o \
	linalg.o \
	jac_rhs.o \
	backward_euler.o \
	nse_solver.o \
	bader_deuflhard.o \
	solver.o \
	decays.o

libsolver.a: $(OBJS) c_fortran_dgssv.o $(patsubst %.o,libsolver.a(%.o),$(OBJS))
	ranlib libsolver.a

c_fortran_dgssv.o:
	ln -sf ./SuperLU_5.0/FORTRAN/c_fortran_dgssv.o c_fortran_dgssv.o

%.o %.mod: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

clean:
	$(RM) *.o *.mod libsolver.a
	$(RM) c_fortran_dgssv.o

# declare module dependencies
include solver.dep
