import matplotlib.pyplot as _plt
import numpy as _np

def draw_network(n,z,lw,c):
    # plot outline of the network, starting with horizontal lines
    zuq = _np.unique(z)
    for zz in zuq:
        ns = _np.sort(n[_np.where(z == zz)])
        seqs = _np.split(ns, _np.where(_np.diff(ns) != 1)[0]+1)
        for seq in seqs:
            y = zz; x1 = _np.min(seq); x2 = _np.max(seq)
            _plt.plot([x1-.5,x2+.5],[y+.5,y+.5],c=c,lw=lw)
            _plt.plot([x1-.5,x2+.5],[y-.5,y-.5],c=c,lw=lw)

    # now the vertical lines
    nuq = _np.unique(n)
    for nn in nuq:
        zs = _np.sort(z[_np.where(n == nn)])
        seqs = _np.split(zs, _np.where(_np.diff(zs) != 1)[0]+1)
        for seq in seqs:
            x = nn; y1 = _np.min(seq); y2 = _np.max(seq)
            _plt.plot([x+.5,x+.5],[y1-.5,y2+.5],c=c,lw=lw)
            _plt.plot([x-.5,x-.5],[y1-.5,y2+.5],c=c,lw=lw)

def draw_magic(lw=.5):
    """draw the neutron and proton magic numbers"""
    magic = [8, 20, 28, 50, 82, 126]
    for m in magic:
        _plt.axvline(m-.5,ls=':',c='grey',lw=lw)
        _plt.axvline(m+.5,ls=':',c='grey',lw=lw)
        _plt.axhline(m-.5,ls=':',c='grey',lw=lw)
        _plt.axhline(m+.5,ls=':',c='grey',lw=lw)


