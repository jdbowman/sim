"""
Generate a SIN (orthographic) projected FITS images from a HEALPix image.
The HEALPix image is assumed to be in eliptical coordinate system. No
coordinate rotation is implemented at the moment.

"""
import argparse
import numpy as np
import healpy as hp
from astropy import wcs
from astropy.io import fits
from multiprocessing import Pool


class _hpm:
    """
    A global object to shared an input HEALPix map among processes to reduce
    the memory comsumption.

    """
    map=None


def _hpm2sin_base(args):
    """
    Base calculation for hpm2sin. Implemented to be called with
    multiprocessing.Pool.map

    """
    fitsfile, ra, dec, dim, res, multiplier, hdr = args
    print 'extracting' + fitsfile

    # Create a new WCS object. The number of axes must be set from the start
    w = wcs.WCS(naxis=2)

    # Set up a SIN projection
    w.wcs.crpix = [dim / 2, dim / 2]
    w.wcs.cdelt = [-res, res]
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

    # Now, write out the WCS object as a FITS header, adding additional
    # fits keyword as applied
    header = w.to_header()
    if hdr is not None:
        for key, value in hdr.iteritems():
            header[key] = value

    # Some pixel coordinates of interest.
    x, y = np.mgrid[0:dim, 0:dim]

    # Convert pixel coordinates to world coordinates
    ra, dec = w.wcs_pix2world(x, y, 0)
    valid_pix = np.isfinite(ra)
    ra *= np.pi / 180
    dec = np.pi * (dec + 90) / 180

    # Get the pixel value from the HEALPix image
    proj_map = np.zeros((dim, dim))
    proj_map[valid_pix] = multiplier * hp.get_interp_val(_hpm.map, dec[valid_pix], ra[valid_pix])

    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file. data is the image array. Axes in 2D numpy
    # array are ordre slow then fast, opposite to fits ordering, so we have to
    # transpose our image array
    hdu = fits.PrimaryHDU(data=proj_map.T, header=header)
    hdu.writeto(fitsfile, clobber=True)


def hpm2sin(hpmfile, fitsfile, ra, dec, dim=7480, res=0.015322941176470588,
            multiplier=1, nthreads=1, hdr=None):
    """
    Generate a SIN (orthographic) projected FITS images from a HEALPix image.
    The HEALPix image is assumed to be in eliptical coordinate system. No
    coordinate rotation is implemented at the moment.

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
    nthreads: integer, optional
        number of processes to be spawned
    hrd: dictionary, optional
        header dictionary to be added or modify to the output fits image.
        format: {'fits keyword' : (value, comment), ...}

    Note
    ----
    This program can project a single HEALPix image to multiple SIN projected
    FITS images if fitsfile, ra, dec and multiplier arguments are given as
    array-like objects. nthreads keyword can be used to parallelize multile
    projection. The default combination of dim and res give a half-sky SIN
    image with ~0.9" resolution, suitable as inputs for MWA simulation in MAPS.

    """
    print "Reading {:s} to a shared memory block"
    _hpm.map = hp.read_map(hpmfile)
    if (isinstance(fitsfile, (np.ndarray, list, tuple)) and
       isinstance(ra, (np.ndarray, list, tuple)) and
       isinstance(dec, (np.ndarray, list, tuple)) and
       isinstance(multiplier, (np.ndarray, list, tuple))):
        args = [(f, r, d, dim, res, m)
                for f, r, d, m in zip(fitsfile, ra, dec, multiplier, hdr)]
    else:
        args = [(fitsfile, ra, dec, dim, res, multiplier, hdr)]
    workers = Pool(nthreads)
    workers.map(_hpm2sin_base, args)


# Command-line paarsing
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