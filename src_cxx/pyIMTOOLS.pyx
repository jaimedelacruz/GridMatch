"""
CYTHON interface for C++ math tools.
Author: J. de la Cruz Rodriguez (ISP-SU, 2019)
"""
import cython
#from cython.parallel import parallel, prange
#import numpy as np
cimport numpy as np
from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray, zeros, ones, outer, int32, arange, nan_to_num, nanmean, float64, where, isfinite, polyfit
from libcpp cimport bool
import sys
from scipy.interpolate import interp2d

__author__="Jaime de la Cruz Rodriguez (ISP-SU 2019)"
__status__="Developing"
__email__="jaime@astro.su.se"


cdef extern from "libgrid.hpp":

    cdef void stretch_matrix_double(int  ny, int  nx, int  npy, int  npx,
			   const double* const  gr, double* const  out_y, double* const out_x, int nthreads)
    
    cdef void stretch_double(int  ny, int  nx, const double* const im, int npy, int npx,
		             const double* const  gr, double* const  res, int nthreads)
    
    cdef void dsgridnest(int  ny, int  nx,  const double* const im1, const double* const  im2,
		         int  nTile,  const int* const  tile,  const int* const clip,
		         double* const res, int nthreads)
    
cdef extern from "imtools.hpp":
    cdef void nearest2D(const int ny, const int  nx, double*  y, double*  x, double*  d_in,
                        const int  ny1, const int  nx1, double* yy_in, double*  xx_in,
                        double*  res_in, const int nthreads, double missing);
    cdef void  bilint2D(const int ny, const int  nx, double*  y, double*  x, double*  d_in,
                        const int  ny1, const int  nx1, double* yy_in, double*  xx_in,
                        double*  res_in, const int nthreads, double missing);
    cdef void  bilint_fast2D(const int ny, const int  nx, double*  d_in,
                        const int  ny1, const int  nx1, double* yy_in, double*  xx_in,
                             double*  res_in, const int nthreads);


def interpolate2D(ar[double,ndim=1] y, ar[double,ndim=1] x, ar[double,ndim=2] im, ar[double,ndim=2] yy, ar[double,ndim=2] xx, int nthreads = 4, nearest=False, missing = None):
    cdef int ny = im.shape[0]
    cdef int nx = im.shape[1]
    
    cdef int ny1 = yy.shape[0]
    cdef int nx1 = yy.shape[1]

    cdef ar[double,ndim=2] res = empty((ny1,nx1), dtype='float64', order='c')

    cdef double miss = 0.0
    
    if(missing is None):
        miss = nanmean(im)
    else:
        miss = <double>missing
    
    if(nearest):
        nearest2D(ny,nx,<double*>y.data, <double*>x.data, <double*>im.data, ny1, nx1, <double*>yy.data, <double*>xx.data, <double*>res.data, <int>nthreads, <double>miss)
    else:
        bilint2D(ny,nx,<double*>y.data, <double*>x.data, <double*>im.data, ny1, nx1, <double*>yy.data, <double*>xx.data, <double*>res.data, <int>nthreads, <double>miss)
        
    return res

def interpolate_fast2D(ar[double,ndim=2] im, ar[double,ndim=2] yy, ar[double,ndim=2] xx, int nthreads = 4):
    cdef int ny = im.shape[0]
    cdef int nx = im.shape[1]
    
    cdef int ny1 = yy.shape[0]
    cdef int nx1 = yy.shape[1]

    cdef ar[double,ndim=2] res = empty((ny1,nx1), dtype='float64', order='c')

    
    bilint_fast2D(ny,nx, <double*>im.data, ny1, nx1, <double*>yy.data, <double*>xx.data, <double*>res.data, <int>nthreads)
        
    return res



def pyStretch(ar im, ar grid, int nthreads=4):

    cdef int ny = im.shape[0]
    cdef int nx = im.shape[1]
    cdef int pny = grid.shape[0]
    cdef int pnx = grid.shape[1]

    cdef ar[double, ndim=2] res = zeros((ny,nx), dtype='float64')
    cdef ar[double, ndim=2] im1
    cdef ar[double, ndim=3] grid1
    
    if(im.dtype == 'float64' and im.dtype == 'float64'):
        stretch_double(ny, nx, <double*>im.data, pny, pnx, <double*>grid.data, <double*>res.data, <int>nthreads)
    else:
        im1 = zeros((ny,nx), dtype='float64')
        grid1 = zeros((pny,pnx,2), dtype='float64')

        im1[:,:] = im[:,:]
        grid1[:,:,:] = grid[:,:,:]
        stretch_double(ny, nx, <double*>im1.data, pny, pnx, <double*>grid1.data, <double*>res.data, <int>nthreads)

    return(res)


def pyStretchMatrix(ar im, ar grid, int nthreads=4, int relative = 0):

    cdef int ny = im.shape[0]
    cdef int nx = im.shape[1]
    cdef int pny = grid.shape[0]
    cdef int pnx = grid.shape[1]

    cdef ar[double, ndim=3] res = zeros((2,ny,nx), dtype='float64')
    cdef ar[double, ndim=3] grid1

    
    if(im.dtype == 'float64' and im.dtype == 'float64'):
        stretch_matrix_double(ny, nx, pny, pnx, <double*>grid.data, <double*>&res[1,0,0], <double*>&res[0,0,0], <int>nthreads)

    else:
        grid1 = zeros((pny,pnx,2), dtype='float64')
        grid1[:,:,:] = grid[:,:,:]
        
        stretch_matrix_double(ny, nx, pny, pnx, <double*>grid1.data, <double*>&res[1,0,0], <double*>&res[0,0,0], <int>nthreads)

    cdef int ii=0
    cdef int jj=0

    # make the correction relative instead of absolute
    if(relative > 0):
        for jj in range(ny):
            for ii in range(nx):
                res[0,jj,ii] -= ii
                res[1,jj,ii] -= jj
    
    return(res)



cpdef congrid2(ar[double,ndim=2] im, int ny, int nx, int nthreads=4):
    cdef int ny1 = im.shape[0]
    cdef int nx1 = im.shape[1]

    #cdef ar[double,ndim=1] xx1 = (arange(nx1, dtype='float64')/(nx1-1.0))
    #cdef ar[double,ndim=1] yy1 = (arange(ny1, dtype='float64')/(ny1-1.0))
    
    cdef ar[double,ndim=2] xx = outer(ones(ny, dtype='float64'), arange(nx,dtype='float64'))* (nx1-1.0) / (nx-1.0)
    cdef ar[double,ndim=2] yy = outer(arange(ny, dtype='float64'), ones(nx,dtype='float64'))* (ny1-1.0) / (ny-1.0)

    cdef ar[double, ndim=2] res = zeros((ny, nx), dtype='float64')
    
    bilint_fast2D(ny1, nx1, <double*>im.data, ny, nx, <double*>yy.data, <double*>xx.data, <double*>res.data, nthreads)


    return res
    
    
def pyGridmatch(ar im1_in, ar im2_in, vgin, clipsin, int nthreads=4):
    """
    pyGridmatch performs a nested search of the distorsion grid that maps im2 into im1.
    
    Arguments:
       im1: Reference image
       im2: Image to map
        vg: numpy array/tuple/list (int) containing the number of tiles to split the image
     clips: numpy array/tuple/list (int) containing the maximum shift allowed within each nested search

    Ported from IDL (dsgridnest.pro) by J. de la Cruz Rodriguez (ISP-SU 2019)

    Bugs found in dsgridnest.pro:
       1) gx and gy are defined as floats but the DLM expects int32
       2) inside ana_gridmatch the images preserve their precision, but all operations should be done as doubles

    """
    cdef int ny = im1_in.shape[0]
    cdef int nx = im1_in.shape[1]

    cdef ar[double,ndim=3] zer = zeros((1,1,2),dtype='float64')
    cdef ar[int,ndim=1] vg     = int32(vgin)
    cdef ar[int,ndim=1] clips  = int32(clipsin)
    cdef int nest = vg.size
    
    if(nest != clips.size):
        print('[error] pyGridmatch: the clips and tile arrays have different dimensions, exiting')
        return zer
    
    
    if((nx != im2_in.shape[1]) or (ny != im2_in.shape[0])):
        print('[error] pyGridmatch: images have different dimensions, exiting')
        return zer
    
    cdef ar[double,ndim=2] im1 = zeros((ny, nx), dtype='float64', order='c')
    cdef ar[double,ndim=2] im2 = zeros((ny, nx), dtype='float64', order='c')

    im1[:,:] = im1_in[:,:]
    im2[:,:] = im2_in[:,:]

    cdef int n = vgin[-1]
    cdef int ngw = int((2*nx)/float(n));
    cdef int nw  = int(1.25 * ngw);

    cdef int ngx =0
    cdef int ngy = 0
    
    if(nx > ny):
        ngx = n
        ngy = int(float(n*ny)/nx +0.5)
    else:
        ngy = n
        ngx = int(float(n*nx)/ny +0.5)

    cdef ar[double,ndim=3] displ = zeros((ngy,ngx,2), dtype='float64')
    
    dsgridnest(ny, nx, <double*>im1.data, <double*>im2.data, nest, <int*>vg.data, <int*>clips.data, <double*>displ.data, nthreads)

    
    
    return displ


# ***********************************************************************

def destretch_tseries(ar cub, double platescale, int tstep, grids=[1], clips=[1], int nthreads=2):

    if(cub.ndim != 3):
        print("[error] destretch_tseries: you must provide a 3D numpy cube, exiting")
        return 0.0
    
    cdef int nt = cub.shape[0]
    cdef int ny = cub.shape[1]
    cdef int nx = cub.shape[2]

    cdef ar[int,ndim=1] Grids = int32(grids)
    cdef ar[double,ndim=1] Clips = int32(clips)
    cdef ar[double,ndim=2] refim = cub[0]

    cdef int maxg      = Grids.max()
    cdef double maxstr = 2*platescale

    cdef int maxgx = maxg
    cdef int maxgy = round(float(ny)/nx*maxgx)
    if(nx <= ny):
        maxgy = maxg
        maxgx = round(float(nx)/ny*maxgy)

    cdef ar[double,ndim=4] delta = zeros((nt,maxgy,maxgx,2), dtype='float64')
    
    cdef int ii = 0
    cdef ar[double,ndim=2] im = zeros((ny,nx), dtype='float64')
    cdef ar[double,ndim=3] dq = zeros((maxgy,maxgx,2), dtype='float64')
    
    for ii in range(1,nt):
        im[:,:] = cub[ii]
        dq[:] = pyGridmatch(refim, im, Grids, Clips, nthreads=<int>nthreads)
        idx = where(dq >  maxstr)
        dq[idx] = 0.0

        delta[ii] = dq
        refim[:] = im

    idx = where(isfinite(delta))
    delta[idx] = 0.0

    # Detrend
    delta = _destretch_gridprep(delta, tstep)

    return delta

# ***********************************************************************

cpdef _smooth(ar[double,ndim=1] var, int win):

    cdef int j0 = 0
    cdef int j1 = 0
    cdef int jj = 0
    cdef int ii = 0
    cdef int N = var.size
    cdef int win2 = win // 2
    
    cdef ar[double,ndim=1] res = zeros(N, dtype='float64')
    cdef double suma = 0.0
    
    for ii in range(N):

        j0 = max(0,ii-win2)
        j1 = min(ii+win2,N)

        suma = 0.0

        for jj in range(j0,j1):
            suma += var[jj]

        res[ii] = suma/max((j1-j0),1)
        
    return res
    
    
# ***********************************************************************

def _destretch_gridprep(delta, int tstep=15):
    cdef int nt = delta.shape[0]
    cdef int ny = delta.shape[1]
    cdef int nx = delta.shape[2]
    cdef int no = delta.shape[3]

    if(no != 2):
        print('[error] _destretch_gridprep: the delta array has non-standard dimensions, exiting')
        return delta


    if(tstep > nt):
        tstep = nt*1
    if(tstep < 3):
        tstep = 3

    cdef ar[double,ndim=2] delx = ascontiguousarray(delta[:,:,:,0])
    cdef ar[double,ndim=2] dely = ascontiguousarray(delta[:,:,:,1])

    cdef int dsz = delx.size
    cdef ar[double,ndim=1] xvec = arange(nt, dtype='float64')
    cdef ar[double,ndim=1] yvec = zeros(nt, dtype='float64')

    cdef int ngx = nx
    cdef int ngy = ny

    cdef int ii = 0
    cdef int jj = 0

    cdef ar[double,ndim=1] xq = zeros((nt), dtype='float64')
    cdef ar[double,ndim=1] yq = zeros((nt), dtype='float64')

    cdef double xq1 = 0.0
    cdef double yq1 = 0.0
    
    cdef ar[double,ndim=1] cf = zeros(2,dtype='float64')
    
    for jj in range(ngy):
        for ii in range(ngx):
            xq[:] = delx[:,jj,ii].cumsum()            
            xq1 = xq[1] # anchor point for first stretched image
            cf[:] = polyfit(xvec,xq,1)
            yvec[:] = xvec*cf[0] + cf[1]
            xq -= yvec
            xq -= _smooth(xq,tstep)
            xq += xq1 - xq[1]
            delta[:,jj,ii,0] = xq

            yq[:] = dely[:,jj,ii].cumsum()            
            yq1 = yq[1] # anchor point for first stretched image
            cf[:] = polyfit(xvec,yq,1)
            yvec[:] = xvec*cf[0] + cf[1]
            yq -= yvec
            yq -= _smooth(yq,tstep)
            yq += yq1 - yq[1]
            delta[:,jj,ii,1] = yq

    delta[0] = 0.0

    return delta

# ***********************************************************************
