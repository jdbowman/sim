import numpy as np
import healpy as hp
import pyfits as pf
import shutil
from subprocess import call
from multiprocessing import Pool


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        msg  -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = msg


def extrapolate_hpm(hpm_in_name, hpm_out_name, nside_out, fitrange,
                    hpm_in_lmax=None, hpm_in_beam_fwhm=None):
    """Extrapolate a healpix map to higher nside by filling it with Gaussian
    noise. Specifically, the program linearly fits the angular power spectrum of
    the input map and exptrapolates its power spectrum from the fit to desired
    multipole number. It then suinmapracts the extropolated power spectrum from the
    original power spectrum. The residual is then converted into a healpix map
    and added back to the input healpix map interpolated to desired nside.

    Parameters
    ----------
    hpm_in_name: str
        Name (and path) of the input healpix map
    hpm_out_name: str
        Name (and path) of the extrapolated healpix map
    extraptol: int
        Multiple number to extrapolate to
    fitrange: tuple of two int
         Range of multipole numbe to fit the power spectrum,
         e.g. fitrange=(min l to fit, max l to fit)
    hpm_in_lmax: int, optional
        Maxmimum multipole number (l) to calculate cl and alm from the map
    hpm_in_beam_fwhm: float, optional
        FWHM of the Gaussian beam of the input map in radian.
        Some full-sky maps has already been smoothed with Gaussian beam with
        a certain FWHM. Thus, the power spectrum of the map is not the "true"
        power spectrum of the sky.
        - NOT IMPLEMENTED YET -

    """
    # Read the input map
    hpm = hp.read_map(hpm_in_name)
    # Calculate the map power spectrum
    hpm_cl, hpm_alm = hp.anafast(hpm, lmax=hpm_in_lmax, alm=True)

    # NOT IMPLEMENTED YET
    # Calculate power spectrum of a Gaussian beam
    # if hpm_in_beam_fwhm is not None:
    #     hpm_bl = hp.sphtfunc.gauss_beam(hpm_in_beam_fwhm, lmax=hpm_in_lmax)
    #     hpm_cl_corrected = hpm_cl / (hpm_bl * hpm_bl)
    # NOTE IMPLEMENTED YET

    # Fitting a power law to Cl
    # We should do fitting on corrected power spectrum, but this will require
    # that we deconvolve the map before adding extrapolated structure. At this
    # point, I do not know how to do that, so we will stick with the uncorrected
    # power spectrum and map.
    l = np.arange(hpm_cl.size)
    try:
        hpm_cl_fit = np.poly1d(np.polyfit(np.log10(l[fitrange[0]:fitrange[1]]),
                               np.log10(hpm_cl[fitrange[0]:fitrange[1]]), 1))
    except TypeError:
        print 'fitrange should be a tuple of two integer for min and max of \
        fitting range'

    # Put together our extrapolated power spectrum
    original_cl = np.concatenate((hpm_cl,
                                 np.zeros(3 * nside_out - 1 - hpm_cl.size)))
    extrap_cl = np.copy(original_cl)
    extrap_cl[fitrange[1]:] =\
        10 ** (hpm_cl_fit(np.log10(np.arange(fitrange[1], extrap_cl.size))))
    diff_cl = -1 * (original_cl - extrap_cl)
    diff_cl[diff_cl < 0] = 0.  # else we have a problem

    # Convert the different in power spectrum to a Gaussian random field
    hpm_diff = hp.synfast(diff_cl, nside_out)

    # Interpolate the input hpm map to higher NSide
    # healpy has a ud_grade function to interpolate map to higher side, but it
    # does not allow clipping of power spectrum before interpolation of the map.
    hpm_original = hp.alm2map(hpm_alm, nside_out)

    # Our diffuse model is the original hpm map interpolated to higher nside,
    # plus the small structure map.
    hpm_extrap = hpm_original + hpm_diff
    hp.write_map(hpm_out_name, hpm_extrap, dtype=np.float64, fits_IDL=False,
                 coord='G')


def freq_interp_fits(infits, infreq, outfreqs, outfits=None, beta=None,
                     gamma=None, T2Jy=True, nthreads=1):
    """
    Interpolate a FITS image to other frequencies usig a spectral index map
    (beta) and a spectral index curverture map (gamma)

    Parameters
    ----------
    infits: string
        A FITS image
    infreq: float
        Frequency of the input image [MHz]
    outfreqs: float, or numpy array of float
        Output frequencies [MHz]
    outprefix: string, or numpy array of string
    beta: string
        Path and name od a spectral index map associated with the inpu image.
        Assume beta = -2.5 on all pixels if None
    gamma: string
        Path and name of a spectral index curverture map associated with the
        input image. Assume curverture = 0 on all pixels if None

     """
    inmap, inmap_header = pf.getdata(infits, header=True)
    basefreq = infreq
    # print 'Read base map from {0:s}'.format(infits)
    # print 'Frequency of inmap: {0:f} MHz'.format(infreq)
    if isinstance(outfreqs, float):
        freqs = [outfreqs]
    elif isinstance(outfreqs, (tuple, list, np.ndarray)):
        freqs = outfreqs
    else:
        raise Exception("Check outfreqs format")
    print 'Output frequencie(s):'
    print freqs
    if isinstance(beta, (int, float)):
        beta = np.ones(inmap.shape) * beta
    elif isinstance(beta, str):
        print 'read beta map from {0:s}'.format(beta)
        beta = pf.getdata(beta)
    elif beta is None:
        beta = np.ones(inmap.shape) * -2.5
        print 'beta is asuume to be -2.5'
    else:
        raise Exception('beta must be a number, a fits image or None for -2.5')
    if isinstance(gamma, (int, float)):
        gamma = np.ones(inmap.shape) * gamma
    elif isinstance(gamma, str):
        print 'read gamma map from {0:s}'.format(gamma)
        gamma = pf.getdata(gamma)
    elif gamma is None:
        gamma = np.zeros(inmap.shape)
        print 'gamma is asuume to be 0'
    else:
        raise Exception('gamma must be a number, a fits image or None for 0')
    if outfits is None:
        outfits = ['{0:s}_{1:.3f}.fits'.format(infits.rsplit('.', 1)[0], freqs[i])
                   for i in range(len(freqs))]
    elif isinstance(outfits, str):
        outfits = ['{0:s}_{1:.3f}.fits'.format(outfits, freqs[i])
                   for i in range(len(freqs))]
    elif isinstance(outfits, (tuple, list, np.ndarray)):
        pass
    else:
        raise Exception('outfits must be a string, an arrylike of string  or None')
    if T2Jy:
        multiplier = [2 * 1.3806488e-23 * 1.e26 * (f * 1.e6) ** 2
                      / (2.99792458e8 ** 2) for f in freqs]
        inmap_header['BUNIT'] = 'Jy/sr'
    else:
        multiplier = np.ones_like(freqs)
    args = (freqs, multiplier, outfits, inmap, inmap_header, beta, gamma, basefreq)
    p = Pool(nthreads)
    p.map(interp, args)


def interp(args):
    f, m, name, inmap, inmap_header, beta, gamma, basefreq = args
    print 'Scaling base map to {0:.3f}MHz and save output to {1:s}'\
        .format(f, name)
    T = m * inmap * np.exp(beta * np.log(f / basefreq)
                           + gamma * (np.log(f / basefreq)) ** 2)
    pf.writeto(name, T, header=inmap_header, clobber=True)


def freq_interp_hpm(inmap, infreq, outfreq, spectral_index=(), curverture=()):
    inmap = hp.read_map(inmap)
    basefreq = infreq
    print 'Read base map from {0:s}'.format(inmap)
    print 'Frequency of inmap: {0:f} MHz'.format(infreq)
    if isinstance(outfreq, float):
        nu = np.array([outfreq])
    elif isinstance(outfreq, np.ndarray):
        nu = outfreq
    else:
        raise Exception("Check outfreq format")
    print 'Output frequencie(s):'
    print nu
    if not spectral_index:
        beta = np.ones(inmap.size) * -2.5
        print 'beta is asuume to be -2.5'
    elif isinstance(spectral_index, str):
        if spectral_index.rsplit('.')[-1] == 'fits':
            beta = hp.read_map(spectral_index)
            print 'read beta map from {0:s}'.format(beta)
    if not curverture:
        gamma = np.zeros(inmap.size)
        print 'gamma is asuume to be 0'
    elif isinstance(curverture, str):
        if curverture.rsplit('.')[-1] == 'fits':
            gamma = hp.read_map(curverture)
            print 'read gamma map from {0:s}'.format(gamma)
    outname = ['{0:s}_{1:.3f}MHz.fits'.format(inmap.rsplit('.', 1)[0], nu[i])
               for i in range(nu.size)]
    for f, name in zip(nu, outname):
        print 'Scaling base map to {0:.3f}MHz and save output to {1:s}'\
            .format(f, name)
        T = np.exp(np.log(inmap) + beta * np.log(f / basefreq)
                   + gamma * (np.log(f / basefreq)) ** 2)
        hp.write_map(name, T, coord='G')


def hpm_extract_facet(args):
    """Wrapper for hpm_extract_facet.py

    """
    infile, outfile, ra, dec, size, res, isys = args
    cmd = ['hpm_extract_facet.py', '-s', '{:.6f}_{:.6f}'.format(ra, dec),
           '--size={:.4f}_{:.4f}'.format(size[0], size[1]),
           '--res={:.16f}'.format(res), '--isys={:s}'.format(isys),
           '-i', infile]
    outname = ('{:s}_{:.6f}_{:.6f}.fits'.format(infile[:-5], ra, dec))
    call(cmd)
    shutil.move(outname, outfile)
    print "Moving {0:s} to {1:s}".format(outname, outfile)


def extract_hpm(infile, outfile, ra, dec, size=(114.6156, 114.6156),
                res=0.9193764705882352, isys='ga', nthreads=1):
    """extract healpix image to a SIN projection FITS image.

    Parameters
    ----------
    infile: str
        Name and path of the healpix image
    outfile: str or array-like of str
        Name and path of the output FITS image
    ra: float or array-like of float
        Right ascension of the center of the FITS image [decimal hours]
    dec: float or array-like of float
        Declination of the center of the FITS image [decimal degree]
    size: array_like of float, 2 elements, optional
        Angular size (x, y) of the output FITS image [degree]
        The default value give a "half" sky image with minimum zero padding.
    res: float, optional
        Angular resolution of the output FITS image. [arcminute]
        The default value gives a FITS image fine enough to be sampled in visgen
        up to 200 MHz

    Note: outfile, ra and dec can be array array-like of equal length to extract
    the input healpix map to multiple FITS files.

    """
    if isinstance(ra, (int, float)):
        ra = [ra]
    if isinstance(dec, (int, float)):
        dec = [dec]
    if isinstance(outfile, str):
        outfile = [outfile]
    if not len(ra) == len(dec) == len(outfile):
        raise InputError("length of outfile, ra and dec are not equal")
    args = [[infile, o, r, d, size, res, isys] for o, r, d in zip(outfile, ra, dec)]
    p = Pool(nthreads)
    p.map(hpm_extract_facet, args)


def scale_fits(infile, outfile, factor):
    hdulist = pf.open(infile)
    hdu = hdulist[1]
    hdu.data.field('TEMPERATURE')[:] = hdu.data.field('TEMPERATURE')[:] * factor
    hdulist[1] = hdu
    hdulist.writeto(outfile, clobber=True)
