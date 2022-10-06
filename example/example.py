import GridMatch as GM
import numpy as np
from astropy.io import fits
import tseries
import matplotlib.pyplot as plt; plt.ion()

# **************************************************

def readFits(filename, ext=1):
    io = fits.open(filename, 'readonly')
    print('reading -> {0}'.format(filename))
    dat = np.ascontiguousarray(io[ext].data, dtype='float32')
    io.close()
    return dat

# **************************************************

if __name__ == "__main__":

    # Let's load a burst of images taken at 35 fps, so the observed
    # object is assumed to be static, and the only source of distorsion
    # is the Earth's atmosphere
    cub = readFits('wb_burst.fits')

    # We define a grid of tiles to split the image
    # the array should go from small to large number of tiles.
    # The routine will perform a nested search using those values.
    # The clips array defines the maximum shift allowed for each tile
    tiles = [8,16,32,48,64]
    clips = [15,10,5,3,2]
    nthreads = 8

    # Calculate the distorsion grid and apply corrections to the
    # cube
    cor, cub1 = GM.DSGridNestBurst(cub, tiles, clips, nthreads = nthreads, apply_correction = True)

    # You can also apply the correction to individual images by hand
    img = GM.Stretch(cub[0], cor[0], nthreads=nthreads)

    
    # show the result, comparing the original cube to the restored one
    nt, ny, nx = cub.shape
    cub2 = np.zeros((nt,ny,nx*2+2), dtype='float32')
    cub2[:,:,0:nx] = cub
    cub2[:,:,nx+2::] = cub1

    tseries.tseries(cub2, dpi = 125)
    
    
    # Now let's assume that we have images where the observed target is slowly
    # evolving (like a CRISP timeseries with a cadence of ~20 sec/frame).
    # In this case we need to use a different routine.
    
    cub3 = readFits('wb_tseries.fits')
    cadence = 30.0 # cadence in seconds
    tstep = round(3*60 / cadence)
    platescale = 1.0 / 0.058 # pixels per arcsec
    
    cor2, cub4 = GM.destretch_tseries(cub3, platescale, tiles, clips, tstep, nthreads = nthreads, apply_correction = True)

    # show the result, comparing the original cube to the restored one
    nt, ny, nx = cub3.shape
    cub5 = np.zeros((nt,ny,nx*2+2), dtype='float32')
    cub5[:,:,0:nx] = cub3
    cub5[:,:,nx+2::] = cub4

    tseries.tseries(cub5, dpi = 125)
    
    
    
    # Now imagine that you want to apply other distorsions
    # (like rotation) and you want to accumulate all the corrections
    # in one matrix, so you only need to interpolate your data
    # once. You can generate matrices for the x,y corrections with
    # the size of the image.

    nt, ny, nx = cub.shape
    corMat = GM.StretchMatrix(ny, nx, cor[0], nthreads=nthreads, bare=True)

    f, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
    ax[0].matshow(corMat[0])
    ax[1].matshow(corMat[1])
    ax[0].set_title('x-axis distorsion matrix')
    ax[1].set_title('y-axis distorsion matrix')
    f.set_tight_layout(True)
