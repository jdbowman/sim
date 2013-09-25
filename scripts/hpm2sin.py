"""
Extract a HEALPix image to half-sky SIN (orthographic) projected FITS images.
The HEALPix image is assumed to be in eliptical coordinate system.

"""
import argparse
import numpy as np
import healpy as hp
from astropy import wcs
from astropy.io import fits
# from multiprocessing.pool import ThreadPool
from multiprocessing import pool
from multiprocessing import sharedctypes
import time

def hpm2sin_base(args):
    print 'Start subprocess'
    start_time = time.time()
    hpm_map_ctypes, shape, fitsfile, ra, dec, dim, res, multiply = args
    hpm_map = np.ctypeslib.as_array(hpm_map_ctypes)
    hpm_map.shape = shape

    print 'extracting' + fitsfile

    # Create a new WCS object. The number of axes must be set from the start
    w = wcs.WCS(naxis=2)

    # Set up a SIN projection
    w.wcs.crpix = [dim / 2, dim / 2]
    w.wcs.cdelt = [-res, res]
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

    # Now, write out the WCS object as a FITS header
    header = w.to_header()

    # Some pixel coordinates of interest.
    x, y = np.mgrid[0:dim, 0:dim]

    # Convert pixel coordinates to world coordinates
    ra, dec = w.wcs_pix2world(x, y, 0)
    valid_pix = np.isfinite(ra)
    ra *= np.pi / 180
    dec = np.pi * (dec + 90) / 180

    # Get the pixel value from the HEALPix image
    proj_map = np.zeros((dim, dim))
    proj_map[valid_pix] = multiply * hp.get_interp_val(hpm_map, dec[valid_pix], ra[valid_pix])

    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file. data is the image array. Axes in 2D numpy
    # array are ordre slow then fast, opposite to fits ordering, so we have to
    # transpose our image array
    hdu = fits.PrimaryHDU(data=proj_map.T, header=header)
    print 'Subprocess finish calculation at {:.6f} seconds'.format(time.time() - start_time)
    hdu.writeto(fitsfile, clobber=True)
    print 'Finish writing fits file at {:.6f} seconds'.format(time.time() - start_time)


def hpm2sin(hpmfile, fitsfile, ra, dec, dim=7480, res=0.015322941176470588,
            multiplier=1, nthreads=1):
    """
    Extract a HEALPix image to half-sky SIN (orthographic) projected FITS images.
    The HEALPix image is assumed to be in eliptical coordinate system.

    Parameters
    ----------
    hpmfile: string
        HEALPix image file
    fitsfile: string or array-like
        output fits file
    ra: float or array-like
        right ascension at the center of the fits image. range=[0,360]deg
    dec: float or array-like
        declination at the center of the fits image. range=[90,-90]deg
    dim: integer, optional
        dimension of the fits image in pixels
    res: float, optional
        angular size of the center pixel in the fits image in deg
    multiplier: float or array-like, optional
        pairs of two numbers, multiplicative and additive, to apply to the fits image

    """
    start_time = time.time()
    hpm_map = hp.read_map(hpmfile)
    size = hpm_map.size
    shape = hpm_map.shape
    hpm_map.shape = size
    hpm_map_ctypes = sharedctypes.RawArray('d', hpm_map)
    hpm_map = np.frombuffer(hpm_map_ctypes, dtype=np.float64, count=size)
    hpm_map.shape = shape
    print 'Finish reading file at {:.6f} seconds'.format(time.time() - start_time)
    if (isinstance(fitsfile, (np.ndarray, list, tuple)) and
       isinstance(ra, (np.ndarray, list, tuple)) and
       isinstance(dec, (np.ndarray, list, tuple)) and
       isinstance(multiplier, (np.ndarray, list, tuple))):
        args = [(hpm_map_ctypes, shape, f, r, d, dim, res, m)
                for f, r, d, m in zip(fitsfile, ra, dec, multiplier)]
    else:
        args = [(hpm_map_ctypes, shape, fitsfile, ra, dec, dim, res, multiplier)]
    # tpool = ThreadPool(nthreads)
    p = pool(nthreads)
    print 'Start passing arguments to workers at {:.6f} seconds'.format(time.time() - start_time)
    # tpool.map(hpm2sin_base, args)
    p.map(hpm2sin_base, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hpmfile', type=str,
        help='HEALPix image file')
    parser.add_argument('fitsfile', type=str,
        help='output FITS file')
    parser.add_argument('ra', type=float, nargs='?', default=0,
        help='right ascension at the center of the fits image. range=[0,360]deg')
    parser.add_argument('dec', type=float, nargs='?', default=0,
        help='declination at the center of the fits image. range=[90,-90]deg')
    parser.add_argument('-d', '--dim', '--dimension', type=int, default=7480,
        metavar='NPIX',
        help='dimension of the fits image in pixels. only support square images.')
    parser.add_argument('-r', '--res', '--resolution', type=float,
        default=0.015322941176470588, metavar='RESOLUTION',
        help='angular size of the center pixel in the SIN projected image')
    parser.add_argument('-m', '--multiply', type=float, default=1,
        metavar='MULTIPLICATIVE',
        help='multiplicative factor to the fits image, ie val_pix *= m')
    parser.add_argument('-l', '--list', action='store_true', metavar='TXTFILE',
        help='threat fitsfile argument as a comma-delimited file containing'
        'with three columns of (filename, ra, dec). Mainly use for multithreading')
    parser.add_argument('--nthreads', type=int, default=1, metavar='NTHREAD',
        help='number of processing threads to use when option l is given')
    args = parser.parse_args()
    if args.list:
        args.fitsfile, args.ra, args.dec = np.genfromtxt(args.fitsfile, unpack=True)
    hpm2sin(args.hpmfile, args.fitsfile, args.ra, args.dec, dim=args.dim,
            res=args.res, multiplier=args.multiply, nthreads=args.nthreads)