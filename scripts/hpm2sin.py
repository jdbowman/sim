"""
Extract a HEALPix image to half-sky SIN (orthographic) projected FITS images.
The HEALPix image is assumed to be in eliptical coordinate system.

"""
import multiprocessing
import argparse
import numpy as np
import healpy as hp
from astropy import wcs
from astropy.io import fits


def hpm2sin_base(args):
    hpm_map, fitsfile, ra, dec, dim, res, scale = args

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
    proj_map[valid_pix] = scale[0] * hp.get_interp_val(hpm_map, dec[valid_pix], ra[valid_pix]) + scale[1]

    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file. data is the image array. Axes in 2D numpy
    # array are ordre slow then fast, opposite to fits ordering, so we have to
    # transpose our image array
    hdu = fits.PrimaryHDU(data=proj_map.T, header=header)
    hdu.writeto(fitsfile, clobber=True)


def hpm2sin(hpmfile, fitsfile, ra, dec, dim=7480, res=0.015322941176470588,
            scale=(1, 0), nthreads=1):
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
    scale: 2-member array-like, optional
        pairs of two numbers, multiplicative and additive, to apply to the fits image

    """
    hpm_map = hp.read_map(hpmfile)
    pool = multiprocessing.Pool(nthreads)
    if (isinstance(fitsfile, (np.ndarray, list, tuple)) and
       isinstance(ra, (np.ndarray, list, tuple)) and
       isinstance(dec, (np.ndarray, list, tuple))):
       args = [(hpm_map, f, r, d, dim, res, scale) for f in fitsfile for r in ra for d in dec]
    else:
        args = (hpm_map, fitsfile, ra, dec, dim, res, scale)
    pool.map(hpm2sin_base, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hpmfile', type=str,
                        help='HEALPix image file')
    parser.add_argument('fitsfile', type=str,
                        help='output FITS file')
    parser.add_argument('ra', type=float,
                        help='right ascension at the center of the fits image. range=[0,360]deg')
    parser.add_argument('dec', type=float,
                        help='declination at the center of the fits image. range=[90,-90]deg')
    parser.add_argument('-d', '--dim', '--dimension', type=int, default=7480, metavar='NPIX',
                        help='dimension of the fits image in pixels. only support square images.')
    parser.add_argument('-r', '--res', '--resolution', type=float,
                        default=0.015322941176470588, metavar='RESOLUTION',
                        help='angular size of the center pixel in the SIN projected image')
    parser.add_argument('-s', '--scale', type=float, nargs=2, default='1 0'.split(),
                        metavar=('MULTIPLICATIVE', 'ADDITIVE'),
                        help='scaling factors to be applied to the fits image, e.g. "-s 3.2 -0.5" will multipliy 3.2 and subtract 0.5 to all pixels in the fits image')
    args = parser.parse_args()
    hpm2sin(args.hpmfile, args.fitsfile, args.ra, args.dec, dim=args.dim,
            res=args.res, scale=args.scale)