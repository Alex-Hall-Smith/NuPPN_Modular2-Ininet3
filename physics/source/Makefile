OBJS = physics_knobs.o \
       reaction_info.o \
       jbj16.o \
       nkk04.o \
       nuc_data.o \
       alpha_decays.o \
       screening_graboske.o \
       screening_chugunov.o \
       screening.o \
       fuller.o \
       reaclib.o \
       nse_torch.o \
       nse_swj.o \
       reverse.o \
       netgen.o \
       kadonis.o \
       vital.o \
       rates.o \
       isomers.o \
       network_boundaries.o \
       networksetup.o \
       abundances.o \
       ppn_physics.o \
       evaluate_rates.o \
       other_nuc.o \
       neutrinos.o

libphysics.a: $(OBJS) $(patsubst %.o,libphysics.a(%.o),$(OBJS))
	ranlib libphysics.a

%.o %.mod: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

clean:
	$(RM) *.o *.mod libphysics.a

# declare module dependencies
include physics.dep
