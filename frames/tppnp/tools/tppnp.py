import numpy as np
import collections
from matplotlib.pyplot import *
from sys import stdout
from fortbin import writeentry
from tqdm import tqdm
import numpy as np
import mmap
import h5py
import resampling
from scipy import interpolate as _interpolate
from scipy.io import FortranFile as _FF
from isonames import *
from abundances import abundances
from progenitor import progenitor
import datetime
import sys as _sys
try:
    import se
except:
    print("Warning: se.py not detected, writing SE files will be disabled")

encd = "ascii"

class particle_set():
    '''
    A set of particles

    Examples:
    ---------
    >>> import tppnp as t
    >>> traj = t.trajectory('/path/to/','trajectory.dat')

    '''
    def __init__(self):
        self._use_offsets               = False
        self._loaded_trajectory         = False
        self._loaded_initial_abundances = False
        self._loaded_final_abundances   = False
        self._loaded_progenitor         = False
        self._filetype = "" # either input or output
        pass

    def __getattr__(self,name):
        if type(self._d) == np.recarray:
            if name not in self._d.dtype.names:
                print("%s not found" % name)
            else:
                return self._d[name]
        else:
            return self._d[name]

    def _init_empty_data(self):
        self._d = {}
        self._d["_pnum"] = np.empty(self.ntracers, dtype=np.int32)
        self._d["_nstp"] = np.empty(self.ntracers, dtype=np.int32)
        self._d["_velfin"] = np.empty((self.ntracers,3), dtype=np.float64)
        self._d["_tpeak"] = np.empty(self.ntracers, dtype=np.float64)
        self._d["_rhopeak"] = np.empty(self.ntracers, dtype=np.float64)
        if self._use_offsets:
            self._d["_offset"] = np.empty((self.ntracers), dtype=np.int64)
        else:
            self._d["_time"] = np.empty((self.ntracers,self.maxnstp), dtype=np.float64)
            self._d["_temp"] = np.empty((self.ntracers,self.maxnstp), dtype=np.float64)
            self._d["_dens"] = np.empty((self.ntracers,self.maxnstp), dtype=np.float64)
            self._d["_dm"]   = np.empty(self.ntracers, dtype=np.float64)
            self._d["_posini"] = np.empty((self.ntracers,3), dtype=np.float64)
            self._d["_posfin"] = np.empty((self.ntracers,3), dtype=np.float64)

    def _assert_can_load_input(self):
        if self._filetype == "output":
            print("This particle set already has type 'output', but you are loading input data.")
            print("Try using a new particle_set object")
            return False
        else:
            return True

    def _assert_can_load_output(self):
        if self._filetype == "input":
            print("This particle set already has type 'input', but you are loading output data.")
            print("Try using a new particle_set object")
            return False
        else:
            return True
        
    def _read_fryer_traj(self):
        """read the data from chris's tracer particles"""
        self.dim = 1 # dimensions of simulation
        f = open(self.tfdir + '/' + self.tfname, 'r')

        # count number of particles in file
        self.ntracers = np.int32(0)
        nstp = 0; maxnstp = 0
        for line in f:
            if 'p' in line:
                self.ntracers += 1
                maxnstp = max(maxnstp,nstp)
                nstp = 0
            else:
                nstp += 1

        # initialise data size
        self.maxnstp = max(maxnstp,nstp)
        self._init_empty_data()

        # For fryer trajectories read the data:
        f.seek(0)
        for i in tqdm(range(self.ntracers),desc="reading trajectories"):
            line = f.readline()
            # tracer attributes
            spl          = line.split()
            self._d['_pnum'][i] = int(spl[0].replace('p',''))
            self._d['_nstp'][i], m, self._d['_dm'][i]  = int(spl[1]), float(spl[2]), float(spl[3])
            self._d['_posini'][i,:] = np.array([m,0,0],np.float64)
            self._d['_posfin'][i,:] = np.array([m,0,0],np.float64)

            # trajectory time-evolution data
            lines = []
            for j in range(self._nstp[i]):
                self._d['_time'][i][j], self._d['_temp'][i][j], self._d['_dens'][i][j] = [float(k) for k in f.readline().split()]

    def _read_townsley_traj(self):
        """read the data from dean's tracer particles"""

        self.dim = 2 # dimensions of simulation
        f = h5py.File(self.tfdir + '/' + self.tfname,'r')

        # read a set of particles from dean's hdf5 file
        particle_idxs = range(0,len(f['trackids']),self.stride)

        self.ntracers = len(particle_idxs)
        self.maxnstp = np.max(f["tracklengths"])
        self._init_empty_data()

        for j in tqdm(range(self.ntracers),desc="reading particle set..."):
            i = particle_idxs[j]
            trackstart  = f["trackstarts"][i]
            tracklength = f["tracklengths"][i]

            # attributes
            self._d['_pnum'][j]       = f["trackids"][i]
            self._d['_dm'][j]         = 1.4e-5
            self._d['_nstp'][j]       = tracklength
            self._d['_posini'][j,0:2] = f["initialpositions"][i]
            self._d['_posfin'][j,0:2] = f["finalpositions"][i]
            self._d['_velfin'][j,0:2] = f["finalvelocities"][i]

            # time evolution data
            self._d['_time'][0:self._nstp[j]] = f["trackdata"][ trackstart:trackstart+tracklength , 0 ]
            self._d['_dens'][0:self._nstp[j]] = f["trackdata"][ trackstart:trackstart+tracklength , 1 ]
            self._d['_temp'][0:self._nstp[j]] = f["trackdata"][ trackstart:trackstart+tracklength , 2 ]


    def _read_LEAFS_traj(self):
        """read the data from LEAFS tracer particles. tfdir == snappath and
        tfname == model. This probably needs to be run on a machine with >200
        GB of RAM for 3d runs.
        Each tracer particle will be recorded from 1 time step before
        the level set arrives, and the temperature in that first time step will
        be set to the temperature at the initial time step in the LEAFS
        simulation (usually 5e5 K)
        """

        try:
            import leafs as l
        except ImportError:
            print("leafs module not found")
            return

        self.dim = 3 # dimensions of simulation
        model = self.tfname
        path = self.tfdir
        t = l.leafs_readtracer(model=model, snappath=path)
        rt = t.loadalltracers()

        ntimesteps, all_times = t.getTimes()

        self.ntracers = len(range(0, t.npart, self.stride))
        self.maxnstp = ntimesteps
        self._init_empty_data()

        particle_ids = range(1, t.npart+1, self.stride)

        for j in tqdm(range(self.ntracers),desc="reading particle set"):
            i = particle_ids[j] - 1
            self._d['_pnum'][j] = particle_ids[j]
            self._d['_dm'][j] = t.tmass[i]
            self._d['_posini'][j,:] = np.array(rt.data[0,i,0:3], np.float64)
            self._d['_posfin'][j,:] = np.array(rt.data[-1,i,0:3], np.float64)
            burning = np.where(rt.data[:,i,14] > 0.)
            if np.size(burning) == 0:
                self._d['_nstp'][j] = 2
                self._d['_time'][j,0:self._nstp[j]] = np.array([all_times[0],all_times[-1]],np.float64)
                self._d['_temp'][j,0:self._nstp[j]] = np.array([5.e5,5.e5],np.float64)
                self._d['_dens'][j,0:self._nstp[j]] = np.array([rt.data[0,i,3],rt.data[-1,i,3]],np.float64)
            else:
                # get trajectory from one time step before the deflagration
                # arrives, setting the temperature back to the initial one
                start = max(np.min(burning) - 1,0)
                self._d['_nstp'][j] = len(all_times[start:])
                self._d['_time'][j,0:self._nstp[j]] = all_times[start:]
                self._d['_temp'][j,0:self._nstp[j]] = np.ascontiguousarray(rt.data[start:,i,4])
                self._d['_temp'][j,0] = rt.data[0,i,4] # usually 5e5 K
                self._d['_dens'][j,0:self._nstp[j]] = np.ascontiguousarray(rt.data[start:,i,3])


    def _read_snsph_traj(self):

        # Read the data from SNSPH 
        self.dim = 3 # dimensions of simulation
        hf = h5py.File(self.tfdir + '/' + self.tfname,'r')

        # get stride info
        selection = range(0,len(hf['particle_id']),self.stride)
        mask = np.zeros(len(hf['particle_id']),dtype='bool')
        mask[selection] = True

        # gather general info about simulation
        self.ntracers = np.count_nonzero(mask)
        particle_idxs = hf['particle_id'][mask]
        self.maxnstp = len(hf['time'])

        self._init_empty_data()

        # store the data to temp memory so you dont reread the h5 file everytime
        mass    = hf['particle_mass'][mask]
        pos_ini = hf[      'pos_ini'][mask,:]
        pos_fin = hf[      'pos_fin'][mask,:]
        time    = hf[         'time']
        rho     = hf[          'rho'][mask,:]
        temp    = hf[         'temp'][mask,:]

        for j in tqdm(range(self.ntracers),desc="reading particle set..."):

            # attributes
            self._d['_pnum'][j] = particle_idxs[j]
            self._d['_nstp'][j] = self.maxnstp
            self._d["_dm"][j]   = mass[j]

            # we dont store vel info yet
            self._d['_posini'][j,0:3] = pos_ini[j,:]
            self._d['_posfin'][j,0:3] = pos_fin[j,:]

            # time evolution data
            self._d['_time'][j,:] = time[:]
            self._d['_dens'][j,:] = rho[j,:]
            self._d['_temp'][j,:] = temp[j,:]



    def load_trajectories(self, tfdir, tfname, source, stride=1000, force=False):
        """
        load in particle data from a simulation including metadata and time evolution 
        but not abundances
        todo: memory mapping for large trajectory files
        
        Parameters:
        -----------
        tfdir : string
            path to trajectory file
        tfname : string
            file name
        source : string
            from what type of simulation are the tracer particles? 'fryer', 'townsley', 'LEAFS', 'SNSPH'
            if 'LEAFS', pass tfdir as snappath and tfname as model
        stride : integer
            sparsity of tracer particles to load from full set
        """

        if not self._assert_can_load_input(): return

        if self._loaded_trajectory and not force:
            print("warning: this particle set already has trajectory information loaded")
            print("if you really want to append trajectories, use force")
            return
            

        self.tfdir  = tfdir
        self.tfname = tfname
        self.stride = stride
        self.source = source

        if self.source == 'fryer':
            self._read_fryer_traj()
        elif self.source == 'townsley':
            self._read_townsley_traj()
        elif self.source == 'LEAFS':
            self._read_LEAFS_traj()
        elif self.source == 'SNSPH':
            self._read_snsph_traj()
        else:
            print('ERROR: Select fryer, townsley, LEAFS, or SNSPH source in trajectory')
            return

        self._loaded_trajectory = True
        self._filetype = "input"

        self.make_particle_map() # create map of particles to indices in arrays
        self.sorted = False # don't sort if loading trajectory data


    def make_particle_map(self):
        """ make map of particle number to array indices"""
        self.particle_map = {}
        for i in range(len(self._pnum)):
            self.particle_map[self._pnum[i]] = i

    def load_tppnp_input(self, fn = "tp.lsb", force = False):
        """ load trajectory and initial abundance data from the tppnp binary input
        fn is the name of the tppnp input binary file
        """

        if not self._assert_can_load_input(): return

        if self._loaded_trajectory or self._loaded_initial_abundances:
            if force:
                self.xiso = []
            else:
                print("trajectories and/or initial abundances already loaded for this set")
                print("use force to overwrite")
                return

        with open(fn,'rb') as fb:

            # read header
            # number of particles in file
            self.ntracers = int(np.fromfile(fb,dtype=np.int32,count=1))

            # max/min steps
            self.minnstp = int(np.fromfile(fb,dtype=np.int32,count=1))
            self.maxnstp = int(np.fromfile(fb,dtype=np.int32,count=1))

            # unfortunately cannot do structured read because of the variable time step
            # data format of file
            # TODO instead save offsets of time-evolution data and abundance
            # data, but read everything else
            self._use_offsets = True
            self._init_empty_data()

            for i in tqdm(range(self.ntracers),desc="loading trajectories & abundances"):
                # particle number
                self._d['_pnum'][i] = int(np.fromfile(fb,dtype=np.int32,count=1))

                # number of time steps
                self._d['_nstp'][i] = int(np.fromfile(fb,dtype=np.int32,count=1))

                # particle mass, initial, final positions and time evolution
                if self._use_offsets:
                    self._d["_offset"][i] = fb.tell()
                    fb.seek(fb.tell() + 7*8 + 8*3*self._nstp[i])
                else:
                    self._d['_dm'][i] = np.fromfile(fb,dtype=np.float64,count=1)
                    self._d['_posini'][i] = np.fromfile(fb,dtype=np.float64,count=3)
                    self._d['_posfin'][i] = np.fromfile(fb,dtype=np.float64,count=3)
                    self._d['_time'][i,0:self._nstp[i]] = np.fromfile(fb,dtype=np.float64,count=self._nstp[i])
                    self._d['_temp'][i,0:self._nstp[i]] = np.fromfile(fb,dtype=np.float64,count=self._nstp[i])
                    self._d['_dens'][i,0:self._nstp[i]] = np.fromfile(fb,dtype=np.float64,count=self._nstp[i])

                # pnum again - skip over
                fb.seek(fb.tell() + 4)

                # isotopes & initial abundances: skip metadata for all but first particle; 
                # i.e. assume all particles have the same isotopes
                if i == 0:
                    self.niso = int(np.fromfile(fb,dtype=np.int32,count=1))
                    if self._use_offsets:
                        pass
                    else:
                        self._d['_xiso'] = np.empty((self.ntracers,self.niso),dtype=np.float64)
                   
                    # isotopic information: A, Z, isotope names
                    # this should be the same for all particles, so assuming that we can skip over
                    # this part of the file with seek, and just read it for the first particle
                    self.isos = []
                    # isotopic information: A and Z
                    self.a = np.fromfile(fb,dtype=np.int32,count=self.niso)
                    self.z = np.fromfile(fb,dtype=np.int32,count=self.niso)

                    for j in range(self.niso):
                        strbits = np.fromfile(fb, dtype=np.int8, count=5)
                        iso = ''.join([chr(item) for item in strbits])
                        self.isos.append(iso)
                    self.isos = np.array(self.isos,dtype="|S5")
                else:
                    # skip over niso, a, z (both 4-byte integers) and isos (5-character string)
                    fb.seek(fb.tell() + 4+(4+4+5)*self.niso)

                # isotopic abundances:
                if self._use_offsets:
                    fb.seek(fb.tell() + 8*self.niso)
                else:
                    self._d['_xiso'][i,:] = np.fromfile(fb,dtype=np.float64,count=self.niso)

        # memory-map the file if using offsets
        if self._use_offsets:
            fb = open(fn,'rb')
            self._mapped_file = mmap.mmap( fb.fileno(), length = 0, access=mmap.ACCESS_READ )
            fb.close()

        self.niso = np.int32(self.niso)

        self._loaded_trajectory = True
        self._loaded_initial_abundances = True
        self.sorted = False # don't sort if loading tppnp input

        self.make_particle_map()


    def load_initial_abundances(self, mode, force=False, **kwargs):
        """ mode = 1:
                provide xini for each mass coordinate: requires mass (cell centred),
                yps, a, z or isos; interpolation will be performed
            mode = 2:
                provide xini for each particle: requires yps(nparticles,nspecies), a, z or isos
            mode = 3:
                provide xini for all particles (all equal): requires yps(nspecies), a, z or isos
        """

        if not self._assert_can_load_input(): return

        if not self._loaded_trajectory:
            print("load particle trajectories with load_trajectories before loading initial abundances")
            return

        if self._loaded_initial_abundances:
            if force:
                pass
            else:
                print("initial abundances already loaded, use force to overwrite")
                return

        # get a and z for each isotope from the ppn naming scheme or vice versa
        if "isos" in kwargs.keys():
            isos = np.array( [ i.encode(encd) if not isinstance(i,bytes) else i for i in kwargs["isos"] ], dtype="|S5" )
            a, z = ppn_name_to_a_z(isos)
        else:
            a = kwargs["a"]; z = kwargs["z"]
            isos = a_z_to_ppn_name(a,z)

        # save a, z and isos
        self.a = np.int32(a); self.z = np.int32(z); self.isos = isos; self.niso = np.int32(len(isos))

        self._d['_xiso'] = np.empty((self.ntracers, self.niso), dtype=np.float64)

        # set the initial abundances for each particle xini
        if mode == 1:
            self._d['_xiso'] = self._interp_x_to_1d_lagrangian_grid(kwargs["mass"], kwargs["yps"])
        elif mode == 2:
            self._d['_xiso'][:,:] = kwargs["yps"][:,:]
        elif mode == 3:
            self._d['_xiso'] = np.repeat(kwargs["yps"],self.ntracers).reshape(self.ntracers,self.niso,order="f")

        self._xiso = np.maximum(self._xiso,1.e-99)
        self._loaded_initial_abundances = True

    def load_final_abundances(self, rundir, force=False, sort=False):
        """ load result from a tppnp run, given output dir, using memory mapping
        
        force: overwrite particle set already in memory
        sort: save argsort result so full arrays may be returned in order
        """

        if not self._assert_can_load_output(): return

        if self._loaded_final_abundances:
            if force:
                self._init_empty_data()
            else:
                print("final abundances already loaded, use force to overwrite")
                return

        # get number of completed particles from file
        ndone = np.loadtxt( rundir + "/num_done.dat", dtype=np.int32 )

        with open( rundir + "/abu.dat", "rb" ) as fb:

            # read header information
            nisos = np.fromfile( fb, np.int32, 1 )[0]
            a     = np.int32( np.fromfile( fb, np.float64, nisos ) )
            z     = np.int32( np.fromfile( fb, np.float64, nisos ) )
            isos  = a_z_to_ppn_name(a,z)

            # get offset of remainder of the file
            offset = fb.tell()

            # memory-map the file
            m = mmap.mmap( fb.fileno(), length = 0, access=mmap.ACCESS_READ )

            # structure of file
            datatype=[('_pnum','<i4'),('_dm','<f8'),('_posini','<f8',(3,)),('_posfin','<f8',(3,)),
                    ('_tpeak','<f8'),('_rhopeak','<f8'),('_xiso','<f8',(nisos,)),
                    ('_flag','|S3')]

            # map data
            self._d = np.ndarray( shape=ndone, dtype=datatype, buffer=m, offset=offset)
            self._d = self._d.view(np.recarray)

            # check for incomplete particles
            if any(self._d._flag != b"EOF"):
                print("warning: the entries for some particles were not complete")

        # store z, a and zuq (unique elemental Zs)
        self.z = np.int32(z); self.a = np.int32(a); self.isos = isos; self.zuq = np.unique(self.z)
        self.niso = np.int32(len(self.isos))

        self._loaded_final_abundances = True

        self.make_particle_map()
        if sort:
            self.sorted = True
            self.sort_index = np.argsort(self._pnum)
        else:
            self.sorted = False

    def load_progenitor(self, prog):
        """load the progenitor star; see progenitor class. give instance of progenitor
        class"""

        if self._loaded_progenitor and not force:
            print("progenitor star already loaded for this particle set.")
            print("use force to overwrite")
            return

        self.star = prog
        self._loaded_progenitor = True

    def write_ppn_input(self,pnum):
        """write particle pnum's trajectory and initial abundances to the files trajectory.input
        and initial_abundance.dat for running with ppn (NuGrid single-zone nucleosynthesis frame.
        Check physics/source/abundances.F90 for iabuini options. this should work with iabuini = 2
        """

        if not self._loaded_trajectory:
            print("trajectory information not yet loaded for particle set")
            return

        if not self._loaded_initial_abundances:
            print("initial abundances not yet loaded for particle set")

        # write trajectory.input
        with open("trajectory.input","w") as f:
            f.write("#\n")
            f.write("# " + "".join([t.ljust(23) for t in ['time','T','rho']]) + "\n")
            f.write("#\n")
            f.write("AGEUNIT".ljust(10) + "SEC\n")
            f.write("TUNIT".ljust(10) + "T9K\n")
            f.write("RHOUNIT".ljust(10) + "CGS\n")
            f.write("ID".ljust(10) + "particle {:8d}\n".format(pnum))
            time = self.time(pnum)
            temp = self.temp(pnum)
            dens = self.dens(pnum)
            for i in range(self.nstp(pnum)):
                f.write("{:23.15e}{:23.15e}{:23.15e}\n".format(
                    time[i],
                    temp[i] / 1.e9,
                    dens[i]
                    )
                    )

        # write initial_abundance.dat
        with open("initial_abundance.dat","w") as f:
            xiso = self.xiso(pnum)
            for i in range(self.niso):
                if type(self.isos[i]) in [bytes,np.bytes_]:
                    isoname = self.isos[i].decode("ascii")
                else:
                    isoname = self.isos[i]
                f.write(isoname + "{:23.15e}\n".format(xiso[i]))

        print("ppn input written to:\n    trajectory.input\n    initial_abundance.dat")


    def _get_idx_for_pnum(self,pnum):
        """find index of particle pnum in dataset""" 
        try:
            return self.particle_map[pnum]
        except KeyError:
            print ( "particle %d does not exist in this dataset" % pnum )
            return -1

    def should_i_sort(self,array):
        if self.sorted:
            return array[self.sort_index]
        else:
            return array


    def pnum(self, p=None):
        if p == None:
            return self.should_i_sort(self._pnum)
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self.should_i_sort(self._pnum[i])

    def nstp(self, p=None):
        if p == None:
            return self._nstp
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self.should_i_sort(self._nstp[i])

    def temp(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 7*8 + 8*self._nstp[i]
            return np.ndarray( shape=self._nstp[i], dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._temp)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._temp[i][:self._nstp[i]]

    def dens(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 7*8 + 2*8*self._nstp[i]
            return np.ndarray( shape=self._nstp[i], dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._dens)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._dens[i][:self._nstp[i]]

    def time(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 7*8
            return np.ndarray( shape=self._nstp[i], dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._time)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._time[i][:self._nstp[i]]

    def dm(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i]
            return np.ndarray( shape=1, dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._dm)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._dm[i]

    def xiso(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 7*8 + 3*8*self._nstp[i] + 4+4+(4+4+5)*self.niso
            return np.ndarray( shape=self.niso, dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._xiso)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._xiso[i]

    def posini(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 8
            return np.ndarray( shape=3, dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._posini)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._posini[i]

    def posfin(self, p=None):
        if self._use_offsets:
            i = self._get_idx_for_pnum(p)
            offset = self._offset[i] + 4*8
            return np.ndarray( shape=3, dtype=np.float64, \
                    buffer=self._mapped_file, offset=offset)
        else:
            if p == None:
                return self.should_i_sort(self._posfin)
            else:
                i = self._get_idx_for_pnum(p)
                if i != -1: return self._posfin[i]

    def mencl(self, p=None):
        """return enclosed mass (only for 1d lagrangian particle sets)"""
        if p == None:
            return np.array([a[0] for a in self.posini()], np.float64)
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self._posini[i][0]

    def tpeak(self, p=None):
        if p == None:
            return self.should_i_sort(self._tpeak)
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self._tpeak[i]

    def rhopeak(self, p=None):
        if p == None:
            return self.should_i_sort(self._rhopeak)
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self._rhopeak[i]

    def velfin(self, p=None):
        if p == None:
            return self.should_i_sort(self._velfin)
        else:
            i = self._get_idx_for_pnum(p)
            if i != -1: return self._velfin[i]

    def upsample(self,*args,**kwargs):
        resampling.upsample(self,*args,**kwargs)

    def downsample(self,*args,**kwargs):
        resampling.downsample(self,*args,**kwargs)

    def yields(self, mode = 0, chop_mass = 1e99, el = False):
        """
        mode 0: total mass of isotopes in final particle set
        mode 1: total mass of isotopes in supernova explosion (including
                progenitor abundances).
        el: True if you want elemental yields, False for isotopic
        chop_mass: delimiting mass between processed and unprocessed ejecta. By default
        it is set to be the edge of the particle set domain, but can be user-specified"""

        # yields from supernova-processed material up to chopmass
        mask = (self.mencl() <= chop_mass)
        SNyield = np.dot(self.dm()[mask], self.xiso()[mask])

        if el:
            mel = np.empty( len(self.zuq), dtype = np.float64 )
            for i in range(len(self.zuq)):
                mel[i] = np.sum(SNyield[(self.z==self.zuq[i])])
            SNyield = mel

        if mode == 0:
            return SNyield
        elif mode == 1:
            if not self._loaded_progenitor:
                print("progenitor star is not loaded")
                return

            # get star yield
            if chop_mass == 1e99: chop_mass = np.max(self.mencl())
            staryield = self.star.yields(mass_cut = chop_mass, el = el)

            # add star yields to SN yields
            if el:
                for i in range(len(self.star.zuq)):
                    idx = np.where(self.zuq == self.star.zuq[i])[0]
                    if len(idx) == 0:
                        print("progenitor element %d not found" % self.star.zuq[i])
                    else:
                        idx = int(idx)
                        SNyield[idx] += staryield[i]
            else:
                SNisos = self.isos
                starisos = self.star.isos()
                for i in range(self.star.nspecies):
                    idx = np.where(SNisos == starisos[i])[0]
                    if len(idx) == 0:
                        print("progenitor isotope %s not found in supernova network\n")
                    else:
                        # catch duplicates (for isomers)
                        idx = int(idx[0]) if len(idx) > 1 else int(idx)
                        SNyield[idx] += staryield[i]
            
            return SNyield

    def isoyield(self, iso="PROT ", chop_mass = 1e99):
        """return yield (in solar masses) of given isotope"""
        if not self._loaded_progenitor:
            # progenitor star not loaded. yields will be only for shock-heated
            # particles
            SNyields = self.yields(mode = 0, el = False)
        else:
            # yield includes progenitor envelope
            SNyields = self.yields(mode = 1, chop_mass = chop_mass, el = False)
        iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso  # encode
        idx = np.where(self.isos == iso)[0]
        if len(idx) == 0:
            print("isotope %s not found" % iso)
            return
        idx = int(idx) if len(idx) == 1 else int(idx[0])
        return SNyields[idx]

    def elyield(self, Z = 1, chop_mass = 1e99):
        """return yield (in solar masses) of given element"""
        SNyields = self.yields(mode = 1, chop_mass = chop_mass, el = True)
        idx = np.where(self.zuq == Z)[0]
        if len(idx) == 0:
            print("isotope %s not found" % iso)
            return
        idx = int(idx) if len(idx) == 1 else int(idx[0])
        return SNyields[idx]

    def radioactive_yields(self):
        """return yields of radioactive isotopes and their reference isotopes,
        in that order"""
        rads = ["AL 26","CL 36","CA 41","TI 44","MN 53","NI 56","FE 60","NB 92", \
                "TC 97","TC 98","PD107","SN126","I 129","CS135","SM146","HF182", \
                "PB205"]
        refs = ["AL 27","CL 35","CA 40","TI 48","MN 55","NI 58","FE 56","MO 92", \
                "RU 98","RU 98","PD108","SN124","I 127","CS133","SM144","HF180", \
                "PB204"]

        yrad = np.array([self.isoyield(r) for r in rads])
        yref = np.array([self.isoyield(r) for r in refs])

        return rads, yrad, refs, yref

    def abundance(self, iso="PROT "):
        """ return final abundance of given iso for every particle"""

        iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso  # encode
        idx = np.where(self.isos == iso)[0]
        if len(idx) == 0:
            print("isotope %s not found" % iso)
            return
        idx = int(idx) if len(idx) == 1 else int(idx[0])
        return self.xiso()[:,idx]

    def write_tppnp_input(self,fn='tp'):
        """
        write input for processing with tppnp

        fn : string
            name of file to write
        """

        if not self._loaded_initial_abundances:
            print("no initial abundance data. load using load_initial_abundances before \
                    writing input for tppnp")
            return

        if _sys.byteorder == "little":
            filename = fn + ".lsb" # little endian
        elif _sys.byteorder == "big":
            filename = fn + ".msb" # big endian

        print(" writing file %s" % filename)
        fb = open(filename, "wb")

        # write particle set metadata
        fb.write( np.int32( self.ntracers ) )
        np.min(self.nstp()).tofile(fb); np.max(self.nstp()).tofile(fb)

        # write particle data
        for i in tqdm(range(self.ntracers),desc="writing tppnp input file"):
            p = np.int32( self.pnum()[i] )
            # write particle number and number of timesteps
            fb.write(p); self.nstp(p).tofile(fb)

            # write mass of particle and initial and final positions.
            # N.B. position is an array of 3 8-byte floats that is either (x, y, 0) [2d Eul],
            # (x, y, z) [3d Eul] or (M_enclosed, 0, 0) [1d, Lag], etc
            self.dm(p).tofile(fb); self.posini(p).tofile(fb); self.posfin(p).tofile(fb)

            # write time, temperature and density
            self.time(p).tofile(fb); self.temp(p).tofile(fb); self.dens(p).tofile(fb)

            # write particle number and number of isotopes
            fb.write(p); self.niso.tofile(fb)

            # write a, z, isonames and mass fractions
            self.a.tofile(fb); self.z.tofile(fb); self.isos.tofile(fb)
            self.xiso(p).tofile(fb)

        fb.close()

    def _interp_x_to_1d_lagrangian_grid(self, mass, yps):
        """Interpolate a given set of abundances, (defined by mass co-ordinate,
        isotope name and abundances) onto the computational grid for this set
        of 1d trajectories. input mass coordinates should be cell-centred"""

        nspecies_in = len(yps[1,:])
        iniabunds = np.zeros([self.ntracers, nspecies_in], np.float64)

        if self.source == "fryer":
            # get array of cell-centred mass coordinates for tracer particles
            # WARNING: assumes tracer particles being mapped to are in ascending
            # mass-coordinate order
            traj_masses = self.mencl().copy()
            traj_masses = np.insert(traj_masses,0,0)
            traj_masses = 0.5*(traj_masses[1:] + traj_masses[:-1])
        else:
            print("_interp_x_to_1d_lagrangian_grid only implemented for source: fryer")
            return

        # interpolate each abundance onto the trajectory grid
        for i in tqdm(range(nspecies_in), desc="interpolating abundances to grid..."):
            # set up interpolant
            finterp = _interpolate.interp1d(mass, yps[:,i])

            # get distribution of abundance on trajectory mass grid
            iniabunds[:,i] = finterp(traj_masses)

        return iniabunds

    def particle2h5part(self,modname,iso="O  16"):
        """
        """
        # find index of isotope
        iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso  # encode
        idx = np.where(self.isos == iso)[0]
        if len(idx) == 0:
            print("isotope %s not found" % iso)
            return
        idx = int(idx) if len(idx) == 1 else int(idx[0])

        mf = self.xiso()[:,idx] # mass fraction of this isotope in each particle
        x = np.empty(len(self.pnum()),dtype=np.float32)
        y = np.empty_like(x); z = np.empty_like(x)
        x[:] = self.posfin()[:,0]; y[:] = self.posfin()[:,1]; z[:] = self.posfin()[:,2]

        # write particles to h5part file
        now = datetime.datetime.now()
        tag = now.strftime("%Y-%m-%d-%H:%M")
        fn = iso.decode("ascii").replace(" ","")+"-"+tag+".h5part"
        f = h5py.File(fn,"w")
        g = f.create_group(modname)
        d1 = g.create_dataset("x",data=x[:], dtype="f")
        d2 = g.create_dataset("y",data=y[:], dtype="f")
        d3 = g.create_dataset("z",data=z[:], dtype="f")
        d4 = g.create_dataset(iso.decode("ascii"),data=mf[:], dtype="f")
        f.close()


    def sph2grid_paraview(self,iso="O  16",stride=5,ncell=5,sl=0.,res=128):
        """
        smooth mass fraction of isotope "iso" for particles onto a grid and write a raw binary
        structured grid file for paraview.
        TODO: make this work for 2d data, too
        stride: stride for particles (i.e. if 1, all particles are done)
        ncell: 1/2 dimension of box around each particle to perform smoothing in
        sl: smoothing length (will be calculated automatically if 0.)
        res: resolution of structured grid to smooth to (i.e. 128 means 128x128x128)
        """
        try:
            from smooth import smooth
        except:
            print("please build smooth module first: `python setup.py build_ext --inplace`")

        # find index of isotope
        iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso  # encode
        idx = np.where(self.isos == iso)[0]
        if len(idx) == 0:
            print("isotope %s not found" % iso)
            return
        idx = int(idx) if len(idx) == 1 else int(idx[0])

        mf = self.xiso()[:,idx] # mass fraction of this isotope in each particle
        x = np.empty(len(self.pnum()),dtype=np.float32)
        y = np.empty_like(x); z = np.empty_like(x)
        x[:] = self.posfin()[:,0]; y[:] = self.posfin()[:,1]; z[:] = self.posfin()[:,2]

        # make grid
        xmax = np.max([np.max(x),np.max(y),np.max(z)])
        xmin = np.max([np.min(x),np.min(y),np.min(z)])
        siz = np.max([abs(xmin),abs(xmax)])
        X = np.linspace(xmin,xmax,res,dtype=np.float32)
        Y = np.linspace(xmin,xmax,res,dtype=np.float32)
        Z = np.linspace(xmin,xmax,res,dtype=np.float32)
        YM, XM, ZM = np.meshgrid(X, Y, Z) # note the odd order numpy returns these in
        mfgrid = np.zeros((res,res,res),dtype=np.float32) # result: mass fraction grid

        # "smooth" uses the Monaghan & Lattanzio SPH kernel, which is 0 above 2 smoothing lengths
        # so set smoothing length to be half the space covered by the number of cells chosen, if
        # sl is not user-defined
        if sl == 0.:
            sl = (ncell*.5)*(X[1] - X[0])

        # make that grid, boi:
        smooth(
                x,y,z,mf, \
                XM, YM, ZM, \
                mfgrid, \
                stride, \
                sl, ncell, res
                )

        # dump it to a raw binary file for paraview to get its chops around
        now = datetime.datetime.now()
        tag = now.strftime("%Y-%m-%d-%H:%M")
        fn = iso.decode("ascii").replace(" ","")+"-"+str(res)+"-"+tag+".raw"
        fb = open(fn,"wb")
        mfgrid.tofile(fb)
        fb.close()
        print("written file %s" % fn)


