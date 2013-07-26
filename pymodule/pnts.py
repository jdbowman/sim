import numpy as np
from subprocess import Popen, PIPE, STDOUT
from mpl_toolkits.basemap import Basemap


def srcgen(prefix, nsources, nside='4096', flux_cut='0.1', seed='0',
           verbose=True):
    """Run srcgen.
    INPUT <trpe>: description [default]
    prefix <str>: prefix of the ouput filename
    nsources <str>: number of total sources to generate
    nside <str>: nside of the base healpix map [4096]
    fluxcut <str>: confusion limit flux_cut [0.1 Jy]
    seed <str>: seed for random number generator [0]
    """
    # Check if file exist
    args = ['srcgen', nside, nsources, seed, flux_cut, prefix]
    run = Popen(args, stdout=PIPE, stderr=STDOUT)
    if verbose:
        for line in run.stdout:
            print(">>> " + line.rstrip())
    stdout, stderr = run.communicate()
    with open(prefix + '.srcgenlog', 'w') as f:
        f.write(stdout)


def mrc2oob(infile, outfile=False):
    """Convert Molonglo source catalog into visgen input format
    """
    # Read RA, Dec amd flux
    cat = np.genfromtxt(infile, comments='#',
                        delimiter=(31, 2, 1, 2, 1, 4, 1, 3, 1, 2, 1, 2, 7, 7),
                        usecols=(1, 3, 5, 7, 9, 11, 13), autostrip=True,
                        unpack=True)
    ra = cat[0] + (cat[1] / 60.) + (cat[2] / 60. / 60.)
    dec = cat[3] + ((cat[4] / 60.) + (cat[5] / 60. / 60.)) * np.sign(cat[3])
    flux = cat[6]
    # Output file
    zero = np.zeros(ra.size, dtype=np.int)  # For QUV rows
    outarray = np.transpose(np.array([ra, dec, flux, zero, zero, zero]))
    if outfile is True:
        prefix = infile.rsplit('/', 1)[-1].rsplit('.', 1)[0]
        outfile = '{0}.oob'.format(prefix)
        np.savetxt(outfile, outarray, fmt='%-10.5f%-12.5f%-8.3f%2d%2d%2d')
    elif isinstance(outfile, str):
        np.savetxt(outfile, outarray, fmt='%-10.5f%-12.5f%-8.3f%2d%2d%2d')


def combine_mrc_srcgen(mrcfile, srcgenfile, outfile=None, scale_to=408,
                       mrc_si=0.8, cutoff=1.0):
    """Combine Molonglo source catalog with a catalog generated from srcgen.
    Output: a catalog in visgen format.
    """
    # Read RA, Dec and flux from the two catalogs
    mrc_ra, mrc_dec, mrc_flux \
        = np.genfromtxt(mrcfile, usecols=(0, 1, 2), unpack=True)
    srcgen_dec, srcgen_ra, srcgen_flux, srcgen_si \
        = np.genfromtxt(srcgenfile, delimiter=',', unpack=True)
    # Sources filtering for MRC
    # Remove rows inwhich the sources have flux below the cut off flux
    mrc_filter = np.where(mrc_flux < cutoff)
    mrc_flux = np.delete(mrc_flux, mrc_filter, 0)
    mrc_ra = np.delete(mrc_ra, mrc_filter, 0)
    mrc_dec = np.delete(mrc_dec, mrc_filter, 0)
    # Flux scaling for MRC
    mrc_flux *= (scale_to / 408.0) ** -mrc_si
    # Source filtering for srcgen
    # Cut srcgen sources located in MRC field with flux > cut off
    srcgen_dec_filter = np.where((srcgen_dec > -85.0) & (srcgen_dec < 18.5))
    srcgen_filter = np.where(srcgen_flux[srcgen_dec_filter] > cutoff)
    srcgen_flux = np.delete(srcgen_flux, srcgen_filter, 0)
    srcgen_ra = np.delete(srcgen_ra, srcgen_filter, 0)
    srcgen_dec = np.delete(srcgen_dec, srcgen_filter, 0)
    srcgen_si = np.delete(srcgen_si, srcgen_filter, 0)
    # Flux scaling for srcgen
    srcgen_flux *= (scale_to / 150.) ** -srcgen_si
    # Stack both catalogs
    ra = np.hstack((mrc_ra, srcgen_ra / 15.))
    dec = np.hstack((mrc_dec, srcgen_dec))
    flux = np.hstack((mrc_flux, srcgen_flux))
    # Output file
    zero = np.zeros(ra.size, dtype=np.int)  # For QUV rows
    outarray = np.transpose(np.array([ra, dec, flux, zero, zero, zero]))
    hdr = ('mrcfile: {0}\n'
           'srcgen file: {1}\n'
           'scale_to: {2}\n'
           'mrc_si: {3}\n'
           'mrc cutoff: {4}\n'
           .format(mrcfile, srcgenfile, str(scale_to), str(mrc_si), str(cutoff)))
    if outfile is None:
        prefix = mrcfile.rsplit('/', 1)[-1].rsplit('.', 1)[0]
        suffix = srcgenfile.rsplit('/', 1)[-1].rsplit('.', 1)[0]
        outfile = '{0}_{1}_{2:.3f}.oob'.format(prefix, suffix, scale_to)
        np.savetxt(outfile, outarray, fmt='%-10.5f%-12.5f%-12.5f%2d%2d%2d',
                   header=hdr)
    elif isinstance(outfile, str):
        np.savetxt(outfile, outarray, fmt='%-10.5f%-12.5f%-12.5f%2d%2d%2d',
                   header=hdr)
    else:
        print hdr
        return outarray


def source_counts(sourcefile, format='oob', bins=100, area=4 * np.pi,
                  flux_multiply=1.0, euclidian=True, log=True):
    """Compute cumulative source counts of a given source catalog.
    OPTION:
    flux_multiply - multiply the read out flux by this factor
    euclidian - multiply the source count by flux ** 2.5 if True
    bins - number of equally spacing bins if bins is an integer
            (in log space if log=True). bins=bin_edges if bins is an ndarray
    log - if true, bin_edges is in log scale
    OUTPUT:
    A tuple of (CDF, bin_center)
    """
    # Check input format
    if format == 'oob':
        flux = np.genfromtxt(sourcefile, usecols=(2,), unpack=True)
    elif format == 'srcgen':
        flux = np.genfromtxt(sourcefile, delimiter=',', usecols=(2,), unpack=True)
    flux *= flux_multiply
    # Making histrogram
    # Genenate bins
    if isinstance(bins, int):
        if log:
            bin_edges = np.logspace(np.log10(np.min(flux)), np.log10(np.max(flux)), num=bins)
        else:
            bin_edges = np.linspace(np.min(flux), np.max(flux), num=bins)
    elif isinstance(bins, np.ndarray):
        bin_edges = bins
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    # Make a histogram
    dN, bin_edges = np.histogram(flux, bin_edges)
    # Count sources starting from max flux
    # dS = bin_edges[1:] - bin_edges[:-1]
    N = np.zeros(bin_centers.size)
    N[::-1] = [np.sum(dN[-i::]) for i in range(1, dN.size + 1)]
    # if euclidian:
        # return ((dN / (dS * area)) * bin_centers ** 2.5, bin_centers)
    # else:
        # return (dN / (dS * area), bin_centers)
    return N, bin_centers


def orthfilter(source, ra0, dec0, outfile, format='oob', FOV='halfsky'):
    """Select point sources from a source catalog that will be observable
    within a specficic FOV on a SIN projection sphere.

    Parameters
    ----------
    source: str
        Name (and path) of the input source catalog.
    ra0: float, or ndarray of float
        RA at the center of the projection in radian
    dec0:, float, ndarray of float
        Dec at the center of the projection in radian
    outfile: str, or ndarray of str
        Name of the filtered catalog.
    format: str
        'oob' or 'srcgen', specifying the format of the source catalog
    FOV: 'halfsky' or float
        Field of view is half sky, or specify the size of a field of view
        in radian

    Example
    -------
    $ import pntstools.orthfilter as orthfilter
    $ cat = 'mrc_srcgen_4096_1.785469e8_150MHz.oob'
    $ ra0 = np.linspace(45, 75, 10) * np.pi / 180.
    $ dec0 = np.ones(10) * -26.0 * np.pi / 180.
    $ outfile = np.array(['cat_inFOV_{0:.2f}d_{1:.2f}d.oob'\
                        .format(i * 180. / np.pi, j * 180. / np.pi)\
                        for i, j in zip(ra0,dec0)])
    $ orthfilter(cat, ra0, dec0, outfile)
    """
    # Check ra0, dec0 and outfile format
    if isinstance(ra0, float) and isinstance(dec0, float)\
            and isinstance(outfile, str):
        ra0 = np.array([ra0])
        dec0 = np.array([dec0])
        outfile = np.array([outfile])
    elif isinstance(ra0, np.ndarray) and isinstance(dec0, np.ndarray)\
            and isinstance(outfile, np.ndarray) and outfile.dtype.kind == 'S'\
            and (ra0.size == dec0.size == outfile.size):
        pass
    else:
        raise Exception('ra0 and dec0 must be floating numbers or ndarray of \
                        the same size')

    print 'Reading {0}'.format(source)
    # Read source catalog
    if format == 'oob':
        ra, dec, flux = np.genfromtxt(source, usecols=(0, 1, 2), unpack=True)
        ra *= 15 * np.pi / 180.
        dec *= np.pi / 180.
    elif format == 'srcgen':
        dec, ra, flux, si = np.genfromtxt(source, delimiter=',',
                                          usecols=(0, 1, 2, 3), unpack=True)
        ra *= np.pi / 180.
        dec *= np.pi / 180.
    else:
        raise Exception('source catalog format must be either `srcgen` or \
                        `oob`')
    # Calculate cosine of the distance from the center of the projection,
    # cos(c). Source will not be visible at a given (ra0, dec0) as a center
    # of the projection if if cos(c) is negative.
    for i in range(ra0.size):
        if FOV == 'halfsky':
            print 'Selecting sources on a half sky centered at \
                ({0:.3f}rad, {1:.3f}rad) ...'.format(ra0[i], dec0[i])
            cosc = np.sin(dec0[i]) * np.sin(dec)\
                + np.cos(dec0[i]) * np.cos(dec) * np.cos(ra - ra0[i])
            in_FOV = np.where(cosc > 0.0)
        elif isinstance(FOV, float):
            print 'Selecting sources on a {0:.2f} centered at \
                ({1:.3f}rad, {2:.3f}rad) ...'.format(FOV, ra0[i], dec0[i])
            in_FOV = np.all([ra > ra0[i] - FOV / 2., ra < ra0[i] + FOV / 2.,
                            dec > dec0[i] - FOV / 2., dec < dec0[i] + FOV / 2.],
                            axis=0)
        print '... saving filtered catalog as {0}'.format(outfile[i])
        if format == 'oob':
            ra_out = ra[in_FOV] * 180. / np.pi / 15.
            dec_out = dec[in_FOV] * 180. / np.pi
            flux_out = flux[in_FOV]
            zero = np.zeros(ra_out.size)
            out = np.array([ra_out, dec_out, flux_out, zero, zero, zero]).T
            np.savetxt(outfile[i], out,
                       fmt='%-12.5f%-12.5f%-12.5f%2d%2d%2d')
        elif format == 'srcgen':
            ra_out = ra[in_FOV] * 180. / np.pi
            dec_out = dec[in_FOV] * 180. / np.pi
            flux_out = flux[in_FOV]
            si_out = si[in_FOV]
            np.savetxt(outfile[i],
                       np.array([dec_out, ra_out, flux_out, si_out]).T,
                       delimiter=',', fmt='%.4f%.4f%.6f%.6f')


def plot_pnts(source, ra_0, dec_0, FOV='fullsky', grid_size=15.):
    ra, dec = np.genfromtxt(source, usecols=(0, 1), unpack=True)
    ra *= 15.
    sky = Basemap(projection='ortho', lon_0=-ra_0, lat_0=dec_0, celestial=True)
    if FOV == 'fullsky':
        sky.drawmeridians(np.arange(0., 360., grid_size))
        sky.drawparallels(np.arange(90, -90, -grid_size))
        sky.drawmapboundary(fill_color='White')
        catx, caty = sky(ra, dec)
        sky.scatter(catx, caty, 3, marker='o', color='Black')
    elif isinstance(FOV, float):
        format_label = lambda deg: '{0:2.0f}h{1:2.0f}m'\
            .format(np.floor(deg / 15.), np.remainder(deg, 15.) * 60. / 15.)
        cnreq = np.array([[ra_0 + FOV / 2., dec_0 - FOV / 2.],
                         [ra_0 - FOV / 2., dec_0 + FOV / 2.]])
        cnrxy = np.array(sky(cnreq[:, 0], cnreq[:, 1])).T
        cenxy = np.array(sky(ra_0, dec_0))
        m = Basemap(projection='ortho', lon_0=-ra_0, lat_0=dec_0,
                    celestial=True, llcrnrx=cnrxy[0, 0]-cenxy[0],
                    llcrnry=cnrxy[0, 1]-cenxy[1], urcrnrx=cnrxy[1, 0]-cenxy[0],
                    urcrnry=cnrxy[1, 1]-cenxy[1])
        m.drawmeridians(np.arange(cnreq[0, 0], cnreq[1, 0], -grid_size),
                        labels=[0, 0, 0, 1], fmt=format_label)
        m.drawparallels(np.arange(cnreq[0, 1], cnreq[1, 1], grid_size),
                        labels=[1, 0, 0, 0])
        m.drawmapboundary(fill_color='White')
        catx, caty = m(ra, dec)
        m.scatter(catx, caty, 3, marker='o', color='Black')
