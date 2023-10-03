import numpy as np
from tqdm import tqdm
cimport numpy as np
cimport cython

cdef extern from "math.h":
        double exp(double x)

cdef extern from "math.h":
        double sqrt(double x)

cdef float square(float x):
    return x*x

cdef float cube(float x):
    return x*x*x

cdef int mymax(int a, int b):
    if a > b:
        return a
    else:
        return b

cdef int mymin(int a, int b):
    if a < b:
        return a
    else:
        return b

cdef find_nearest(a,v):
    return np.abs(a-v).argmin()

DTYPE_f32 = np.float32
DTYPE_int32 = np.int32

ctypedef np.int32_t DTYPE_int32_t
ctypedef np.float32_t DTYPE_f32_t

@cython.boundscheck(False)
def smooth(np.ndarray[DTYPE_f32_t,ndim=1] x, np.ndarray[DTYPE_f32_t,ndim=1] y,
        np.ndarray[DTYPE_f32_t,ndim=1] z, np.ndarray[DTYPE_f32_t,ndim=1] f,
        np.ndarray[DTYPE_f32_t,ndim=3] XM, np.ndarray[DTYPE_f32_t,ndim=3] YM,
        np.ndarray[DTYPE_f32_t,ndim=3] ZM, np.ndarray[DTYPE_f32_t,ndim=3] FM,
        int stride, float sl, int ncell, int res):
    """particle set with (x[:], y[:], z[:]) and scalar f[:] to grid (meshgrid)
    with (XM[:,:,:], YM[:,:,:], ZM[:,:,:]), using gaussian with smoothing
    length sl. The scalar is smoothed onto the grid FM"""

    cdef DTYPE_f32_t alfa = 3. / (2.*3.141592653589793)
    cdef DTYPE_f32_t h = sl
    cdef DTYPE_f32_t h2 = h*h
    cdef DTYPE_f32_t two_thirds = 2./3.
    cdef DTYPE_f32_t sixth = 1./6.
    cdef DTYPE_f32_t alf_ov_hn = alfa / (h*h*h)
    cdef DTYPE_f32_t denom = 1. / (2.*3.141592653589793*sl**2)
    cdef DTYPE_f32_t expdenom = 1. / (2.*sl**2)

    cdef DTYPE_int32_t i, j, k, l, ix, iy, iz, xmin, xmax, ymin, ymax, zmin, zmax
    cdef DTYPE_f32_t rad, rad2, s, s2, s3, func

    # loop over all particles
    cdef DTYPE_int32_t length = len(x)

    for l in tqdm(range(0,length,stride),desc="smoothing"):

        # find nearest grid cell to anchor our smoothing around:
        ix = find_nearest(XM[:,0,0],x[l]); iy = find_nearest(YM[0,:,0],y[l]);
        iz = find_nearest(ZM[0,0,:],z[l])

        xmin = mymax(ix-ncell,0); xmax = mymin(ix+ncell,res-1)
        ymin = mymax(iy-ncell,0); ymax = mymin(iy+ncell,res-1)
        zmin = mymax(iz-ncell,0); zmax = mymin(iz+ncell,res-1)
    
        for i in range(xmin,xmax):
            for j in range(ymin,ymax):
                for k in range(zmin,zmax):
                    rad2  = square(XM[i,j,k] - x[l]) + \
                            square(YM[i,j,k] - y[l]) + \
                            square(ZM[i,j,k] - z[l]) 
                    s2 = rad2/h2
                    #func = exp( -rad2 * expdenom )
                    if s2 >= 0 and s2 < 1:
                        s = sqrt(s2)
                        s3 = s2*s
                        func = two_thirds - s2 + 0.5*s3
                        FM[i,j,k] += f[l] * func
                    elif s2 >= 1 and s2 < 4:
                        s = sqrt(s2)
                        func = sixth * cube(2 - s)
                        FM[i,j,k] += f[l] * func
                    else:
                        pass



