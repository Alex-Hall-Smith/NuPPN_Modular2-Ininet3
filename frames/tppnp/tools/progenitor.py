# This class stores info about progenitor models for a supernova particle set.
# Right now they are 1d lagrangian profiles but more could be added.

import numpy as _np
from tqdm import tqdm
from isonames import *

class progenitor():
    def __init__(self, **kwargs):
        """instantiate a progenitor star's abundances/structure, requires
        mass coordinate "m" (nzones), abundances "yps"(nzones,nspecies)
        and either isotope names or A and Z (nspecies).
        Progenitor data is assumed to be ordered in direction of increasing
        radius and to begin at the centre of the star.
        The mass co-ordinate is assumed to be edge-centred"""

        if "file" in kwargs.keys():
            d = np.load(kwargs["file"],mmap_mode='r')
            data = d
        else:
            data = kwargs

        for d in ["m", "yps"]: assert d in data.keys()

        # store input data
        self.m = data["m"]
        self.nzones = len(self.m)
        self.x = data["yps"]
        self.nspecies = _np.shape(self.x)[1]

        if "rho" in data.keys(): self.rho = data["rho"]
        if "t" in data.keys(): self.t = data["t"]
        if "r" in data.keys(): self.r = data["r"]

        # check spatial dimensions match up
        assert _np.shape(self.x)[0] == self.nzones

        # get a and z for each isotope if only isotope names are given
        if "isos" in data.keys():
            print("getting a and z...")
            isos = data["isos"]
            isos = _np.array( [ i.encode(encd) if not isinstance(i,bytes) else i for i in isos ] )
            a, z = ppn_name_to_a_z(isos)
        else:
            for d in ["a", "z"]: assert d in data.keys()
            a = data["a"]; z = data["z"]
        
        # check abundance dimensions match up
        assert len(a) == self.nspecies
        assert len(z) == self.nspecies

        self.a = a; self.z = z

    def dm(self):
        """return Delta(edge-centred mass coordinates), i.e. the mass of each
        zone"""
        # get mass of each zone from edge-centred mass coordinate
        extmass = _np.insert(self.m,0,0)
        return _np.diff(extmass)

    def zuq(self):
        """return unique Z (i.e. elements in set)"""
        return _np.unique(self.z)

    def isos(self):
        """return isotope names"""
        return a_z_to_ppn_name(self.a,self.z)

    def yields(self, mass_cut = 0., el = False):
        """return the abundances above the mass cut. This could be a physical mass cut or the
        dimiting mass between the material processed/not processed by the supernova shock.
        Returns integrated yield of each isotope (Y = sum_i(dm_i*x_ij)) in solar masses.
        If el == True, elemental yields are returned.
        """

        # find index of first cell in progenitor above mass cut
        idx = _np.abs(self.m - mass_cut).argmin()
        if self.m[idx] < mass_cut: idx += 1

        # split first shell precisely at mass cut:
        my_dm = np.insert(self.dm()[idx+1:], 0, self.m[idx]-mass_cut)
        yields = _np.dot(my_dm,self.x[idx:])
        if el:
            mel = np.empty( len(self.zuq), dtype = np.float64 )
            for i in range(len(self.zuq)):
                mel[i] = np.sum(yields[(self.z==self.zuq[i])])
            yields = mel

        return yields

    def abundance(self, iso = "", mass_cut = 0.):
        """Return mass fractions of isotope as a function of mass coordinate.
        If no isotope name is given, all isotopes are returned"""

        # find index of first cell in progenitor above mass cut
        midx = _np.abs(self.m - mass_cut).argmin()
        if self.m[midx] < mass_cut: midx += 1

        if iso == "":
            return self.m[idx:], self.x[idx:]
        else:
            iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso  # encode
            idx = _np.where(self.isos() == iso)[0]
            if len(idx) == 0:
                print("isotope %s not found" % iso)
                return
            idx = int(idx) if len(idx) == 1 else int(idx[0])
            return self.m[midx:], self.x[midx:,idx]

