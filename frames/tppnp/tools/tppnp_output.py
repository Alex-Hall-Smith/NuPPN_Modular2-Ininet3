import numpy as np
import abundances

class tppnp_output():
    '''
    final abundance data for a TPPNP run.

    Parameters:
    -----------
    rundir : string
	path to tppnp run directory. The run directory should contain the
	directories "output" and "output_decay", each containing a .npz file
        with the compressed data

    Examples:
    ---------
    >>> import tppnp as t
    >>> traj = t.trajectory('../tppnp_data/trajectory.dat')
    >>> run = t.tppnp_run(traj,'.')

    Explanation of data quantities:
    -------------------------------
    self.niso  : number of isotopes (int)
    self.nel   : number of elements (int)
    self.pname  : particle names (str, dim [self.nps])
    self.pnum   : particle numbers (int, dim [self.nps])
    self.m      : lagrangian coordinate of tracer particle (float, dim [self.nps])
    self.dm     : mass of tracer particle (float, dim [self.nps])
    self.xiso  : isotopic mass fractions in each tracer; 2D data (float, dim [self.nps, self.nisos])

    '''
    def __init__(self, rundir='.'):
        rundir += '/output/' ; fn = 'abu.dat'
        # load data file
        self._load_abudata(rundir + fn)

        # elemental masses
        self.miso = np.dot(self.dm_all(), self.xiso_all())
        self.mel  = np.array([sum(self.miso[i]) for i in [np.where(self.z==j)[0] for j in self.zuq]])

    def _load_abudata(self, fn):
        '''
        load the final abundance data from the tppnp output file abu.dat
        '''
        f = open(fn, 'r')

        # read number of trajectories
        self.ntracers = 0
        for line in file:
            if 'p' in line:
                self.ntracers += 1

        f.seek(0)

        self.data["pnum"]  = np.empty(ntraj, np.int32)
        self.data["pname"] = np.chararray(ntraj, 8)
        self.data["abund"] = []
        self.data["dm"]    = np.empty(ntraj)
        # read abundance data
        for i in tqdm(range(ntraj)):
            spl = f.readline().split()
            pnum, pname, niso, dm = int(spl[1]), 'p'+str(pnum), int(spl[2]), float(spl[3])
            self.data["pnum"][i] = pnum; self.data["pname"][i] = pname
            # initialise arrays with npecies dimensions
            a, z = np.empty((2,niso), np.int32) ; xiso = np.empty(niso, np.float64)
            for j in range(niso):
                spl = f.readline().split()
                a[j], z[j], xiso[j] = int(spl[0]), int(spl[1]), float(spl[3])

            self.data["abund"].append(abundances(a,z,xiso))

    def _get_idx_for_pnum(self,pnum):
        """find index of particle pnum in dataset""" 
        idx = np.where( self.data["pnum"] == pnum )[0]
        # check that this particle exists
        if len(idx) == 0:
            print ( "particle %d does not exist in this dataset" % pnum )
            return -1
        else:
            return int(idx)

    def dm(self,p):
        i = self._get_idx_for_pnum(p)
        if i != -1: return self.data["dm"][i]

    def pnum(self,p):
        i = self._get_idx_for_pnum(p)
        if i != -1: return self.data["pnum"][i]
        
    def pname(self,p):
        i = self._get_idx_for_pnum(p)
        if i != -1: return self.data["pname"][i]

    def xiso(self,p):
        i = self._get_idx_for_pnum(p)
        if i != -1: return self.data["xiso"][i]

    def dm_all(self):
        return np.array( [ self.dm(p) for p in self.data["pnum"] ] )

    def xiso_all(self):
        return np.array( [ self.xiso(p) for p in self.data["pnum"] ] )

    def plot_elemental_yield(self):
        '''
        plot the ejected mass of each element
        '''
        plot(self.zuq,self.mel,'o-')
        xlabel('Z')
        ylabel('ejected mass / M$_\odot$')




