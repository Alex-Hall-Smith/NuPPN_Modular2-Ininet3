import numpy as np
import inspect
import os
from drawnetwork import *
import matplotlib.pyplot as _plt
from isonames import a_z_to_ppn_name, ppn_name_to_a_z

encd = "ascii"
#encd = "utf-8"

import matplotlib.pyplot as _plt

class abundances():
    def __init__(self,a,z,xiso):
        self._a    = a
        self._z    = z
        self.xiso = xiso
        self.niso = len(a)

    @property
    def a(self):
        return np.int32(self._a)

    @property
    def z(self):
        return np.int32(self._z)

    def isos(self):
        """return isotope names in this abundance set"""
        return a_z_to_ppn_name(self.a,self.z)

    def zuq(self):
        """return unique Z numbers (i.e. the elements)"""
        return np.unique(self.z)

    def els(self):
        """return list of elements in this set of abundances"""
        _, idx = np.unique(self.z, return_index = True)
        self.els = self.isos[idx]
        self.els = np.array([el[:2].capitalize() for el in self.els])

    def n(self):
        """neutron number"""
        return self.a - self.z

    def isoidx(self, a, z):
        '''
        return index of isotope in the self.z, self.z, self.isos, self.misos,
        self.xisos[i,:] arrays
        '''
        wh = np.where((self.a == a) & (self.z == z))[0]
        if len(wh) == 0:
            print('isotope not found')
        else:
            return wh[0]

    def ye(self):
        # neutron excess (nex) and electron fraction (ye)
        nex = sum(( self.a - 2. * self.z) * self.xiso / self.a)
        ye = ( 1. - nex ) * 0.5
        return ye

    def xel(self):
        # calculate elemental abundances
        zuq = self.zuq()
        xel = np.empty( len(zuq), dtype = np.float64 )
        for i in range(len(zuq)):
            xel[i] = np.sum(self.xiso[(self.z==zuq[i])])
        return xel

    def plot_xel(self,**kwargs):
        _plt.plot(self.zuq(), self.xel(), ls='', marker='o', **kwargs)
        _plt.xlabel("$Z$") ; _plt.ylabel("$X$")

    def top(self,n):
        """return top n mass fractions"""
        return sorted(zip(self.xiso,self.isos()))[::-1][:n]

    def plot_chart(self):
        """plot an isotopic chart of abundances"""
        abdir = os.path.dirname(inspect.getfile(inspect.currentframe()))
        data = np.load(abdir + "/stable.npz")
        stablez = data["arr_0"]
        stablen = data["arr_1"]
        zuq = np.unique(self.z); nuq = np.unique(self.n())

        # plot abundances
        x = np.ones((np.max(nuq)+1,np.max(zuq)+1)) * np.nan
        for i in range(len(self.z)):
            x[self.n()[i],self.z[i]] = np.log10(self.xiso[i])

        x = np.ma.masked_where(np.isnan(x),x)
        cm = _plt.pcolormesh(np.arange(int(np.max(nuq)+1))-.5,np.arange(int(np.max(zuq)+1))-.5,x.T,rasterized=True)
        cm.cmap.set_under(color="w")
        cbar = _plt.colorbar()
        cbar.set_label("$\log(X)$")

        draw_network(self.n(),self.z,0.5,"#DCDCDC")
        draw_network(stablen,stablez,0.5,"k")
        draw_magic()

        _plt.xlabel("neutron number")
        _plt.ylabel("proton number")
        _plt.gca().set_aspect(1)
        _plt.draw()


