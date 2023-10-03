# import tppnp module
import tppnp as t

# initialise an instance of particle_set
p = t.particle_set()

# read particle data
p.load_trajectories('input_data/','townsley_etal_2016_fullset.hdf5','townsley')

# load in initial abundance data (same for all particles)
yps = np.array( [0.5,0.5], np.float64 ) ; isos = ['C  12','O  16']
p.load_initial_abundances( 3, yps=yps, isos=isos)

# write input for tppnp to post-process
p.write_tppnp_input()

# load result from tppnp
p.load_final_abundances('.')
