# Loading a set of CCSN tracer particles and initial abundances
# -------------------------------------------------------------

import tppnp as t

p = t.particle_set()

ccsn_dir  = "./input_data"
ccsn_file = "traj15f1_216M1.3aal.dat" # tracer particles from Chris's 1d models
comp_file = "./input_data/M15h.20180320.npz" # composition file

# load tracer particle time evolution from Fryer's CCSN explosions, e.g.:
p.load_trajectories(ccsn_tfdir,ccsn_file,source="fryer")

# load the progenitor abundances and edge-centred mass coordinate from numpy binary file
d = np.load(comp_file)

# get cell-centred mass coordinates
m = np.insert(d['m'],0,0)
m = 0.5 * (m[1:] + m[:-1])

# load initial abundance data into particle set
p.load_initial_abundances(1,mass=m,yps=d['yps'],a=d['a'],z=d['z'])

# write input files for tppnp code
p.write_tppnp_input()

# load result from tppnp once computation is done
p.load_final_abundances('.')

# Changing the rates for the sensitivity study
# --------------------------------------------
# this can now be done with the physics input deck array variables:
#     rate_index(1:10)
#     rate_factor(1:10).
# For example, if the rate we want to modify is number 1001 (check its number
# in networksetup.txt), we can give it a rate factor of 10 by adding the following
# into ppn_physics.input:
#     rate_index(1) = 1001
#     rate_factor(1) = 10.d0
# If detailed balance = .true. is set in the physics input deck as well, the
# reverse rate will be calculated by detailed balance, so you don't need to set
# the reverse rate factor. Be careful, though, if you want to modify a reverse
# rate specifically and you have detailed balance switched on.

