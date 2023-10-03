import numpy as np
from scipy import interpolate
from tqdm import tqdm

def upsample(t, pnums = [], T9lo = 0.3, delta_temp_max = 0.05, \
                 delta_dens_max = 0.05):
        """
        Add time-resolution to a trajectory above a given temperature using
        linear interpolation. Resolution is added based on a maximum allowed
        fractional change in temperature and density

        Parameters:
        -----------
        pname : string
            name of trajectory to upsample; if 'all', upsample all
            trajectories
        T9lo : float
            temperature above which to increase the resolution
        delta_temp_max : float
            max. allowed fractional change in temperature
        delta_dens_max : float
            max. allowed fractional change in density
        """

        if t.source != 'fryer': raise Exception("upsampling only for fryer data")

        tpdo = t.data["pnum"] if pnums == [] else pnums
        
        for i in tqdm(range(len(tpdo))):
            p = tpdo[i]
            time = t.time(p).copy()
            temp = t.temp(p).copy()
            dens = t.dens(p).copy()

            # create interpolants
            temp_interp = interpolate.interp1d(time, temp, kind='linear')
            dens_interp = interpolate.interp1d(time, dens, kind='linear')

            # find locations to insert more steps:
            tempdif = np.abs(np.diff(temp) / temp[:-1])
            densdif = np.abs(np.diff(dens) / dens[:-1])
            nstptemp = tempdif / delta_temp_max
            nstpdens = densdif / delta_dens_max
            nstp = np.maximum(nstptemp, nstpdens)

            # where to insert and how many steps
            avT9 = (temp[1:] + temp[:-1]) / 2. / 1.e9
            tofix = np.where((avT9 >= T9lo) & (nstp > 1))[0]
            nstp = np.array(nstp[tofix], np.int64)
            nstp = nstp + 1 # round up
            
            # insert steps where necessary:
            for i in range(len(tofix)):
                a = time[tofix[i]]; b = time[tofix[i]+1]
                instime = np.linspace(a, b, nstp[i] + 2)[1:-1]
                instemp = temp_interp(instime)
                insdens = dens_interp(instime)
                time = np.insert(time, tofix[i]+1, instime)
                temp = np.insert(temp, tofix[i]+1, instemp)
                dens = np.insert(dens, tofix[i]+1, insdens)
                # now shift indices for inserting by number of steps just
                # inserted, so we stay at the right place
                tofix[i:] = tofix[i:] + nstp[i]

            idx = t._get_idx_for_pnum(p)
            t.data [ "time"   ]  [ idx ]  = time
            t.data [ "temp"   ]  [ idx ]  = temp
            t.data [ "dens"   ]  [ idx ]  = dens
            t.data [ "nstp"   ]  [ idx ]  = len(time)

def downsample(t, pnums = [], T9lo = 0.3, delta_temp_loT = 0.2):
    """
    Reduce time resolution of a trajectory below a given temperature based
    on an ideal fractional change in temperature

    Parameters:
    -----------
    pname : string
        name of trajectory to upsample; if 'all', downsample all
        trajectories
    T9lo : float
        temperature below which to reduce the resolution
    delta_temp_max : float
        max. allowed fractional change in temperature
    """

    if t.source != 'fryer': raise Exception("downsampling only for fryer data")

    tpdo = t.data["pnum"] if pnums == [] else pnums

    for i in tqdm(range(len(tpdo))):
        p = tpdo[i]
        time = t.time(p).copy()
        temp = t.temp(p).copy()
        dens = t.dens(p).copy()

        # downsample the low temperature steps:
        num_to_remove = 1e99 # init
        while num_to_remove > 0:
            avT9 = (temp[1:] + temp[:-1]) / 2. / 1.e9
            tempdif = np.abs(np.diff(temp) / temp[:-1])
            rmvtemp = (tempdif < delta_temp_loT) & (avT9 < T9lo)
            num_to_remove = np.count_nonzero(rmvtemp)
            temp = np.delete(temp, np.where(rmvtemp)[0][::2])
            time = np.delete(time, np.where(rmvtemp)[0][::2])
            dens = np.delete(dens, np.where(rmvtemp)[0][::2])

        idx = t._get_idx_for_pnum(p)
        t.data [ "time"   ]  [ idx ]  = time
        t.data [ "temp"   ]  [ idx ]  = temp
        t.data [ "dens"   ]  [ idx ]  = dens
        t.data [ "nstp"   ]  [ idx ]  = len(time)


