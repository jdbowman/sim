import numpy as np
import pyfits as pf
import pywcs as wcs
import scipy.interpolate as interp


def read_sim_ascii(simfile, sep=',', shape=(), cols=None):
    """
    Read sky coordinate (RA, Dec), frequencies and temperature from 21 cm
    simulation cube.

    Parameters
    ----------
    simfile: string
        A simulation file in ascii table format
    delimiter: string, optional
        Delimiter in the table. [default = ',']
    shape: tuple of intiger, optional
        Shape of the simulation cube. You can assign any valid numpy ndarray
        shape to reshape the read out from the simulation file, else the output
        is a straight array read out from the input file. Note that the program
        only supports the same shape for RA, Dec, frequencies and temperature.
    cols: tuple of integer or None, option
        Columns number to read out RA, Dec, frequencies and temperature
        from the input file.
        e.g. (3,0,1,2) for table formatted to temperature, RA, Dec, frequencies
        None assume sequencial order

    Returns
    -------
    out: ndarray
        RA, Dec, frequencies and temperatures read from the input file

    """
    ra, dec, freq, temp = np.genfromtxt(simfile, delimiter=sep, usecols=cols,
                                        unpack=True)
    ra.shape, dec.shape, freq.shape, temp.shape = shape, shape, shape, shape
    return ra, dec, freq, temp


def read_sim_bin(simfile, dtype=float):
    """
    Read 21 cm simulation cube from a binary file

    Parameters
    ----------
    simfile: string
        A simulation file in ascii table format. Assume "C" order
    dtype: numpy or python type
        type of input binary file

    Returns
    -------
    out: ndarray
        the datacube

    """
    cube = np.fromfile(simfile, dtype=dtype)
    return cube


def read_map(mapfile):
    """Read sky coordinates (RA, Dec) and pixel data of a map (FITS image)
    created by hpm_extract_facet

    Parameters
    ----------
    mapfile: str
        A FITS image of a map created by hpm_extract_facet

    Returns
    -------
    out: ndarray, shape (npoints,)
        RA, Dec and pixel data read from the map
    """
    # Read FITS image
    hdulist = pf.open(mapfile)
    data = hdulist[0].data.squeeze()
    priheader = hdulist[0].header
    # Modify the header because we squeeze the data
    if priheader['NAXIS'] != 2:
        priheader['NAXIS'] = 2
        priheader.pop('NAXIS3')
        priheader.pop('NAXIS4')
    # Read WCS info and convert pixel to sky coordinate
    map_wcs = wcs.WCS(priheader)
    pix_x, pix_y = np.mgrid[0:data.shape[0], 0:data.shape[1]]\
        .reshape((data.ndim, data.size))
    ra, dec = map_wcs.wcs_pix2sky(pix_x, pix_y, 0)
    return ra, dec, data.reshape(data.size, order='F')


def map2sim(mapfile, simfile, outfile, interp_method='nearest'):
    """Read and interpolate a sky map onto 21 cm simulation cube

    Parameters
    ----------
    mapfile: str
        A FITS image od a map created by hpm_extract_facet
    simfile: str
        A simulation cube file
    outfile: str
        Path and name of an output file. outfile is in a simulation cube
        format with interpolated data from the map as an additional column

    interp_method : {'linear', 'nearest', 'cubic'}, optional
        Method of interpolation. One of
        - ``nearest``: return the value at the data point closest to
        the point of interpolation.

        - ``linear``: tesselate the input point set to n-dimensional
        simplices, and interpolate linearly on each simplex.

        - ``cubic`` (1-D): return the value detemined from a cubic spline.
    """
    # Read map and simulation cube
    map_ra, map_dec, map_data = read_map(mapfile)
    sim_ra, sim_dec, sim_freq, sim_data = read_sim_ascii(simfile)
    # Convert RA and Dec of the diffuse map into a cube-centric systemmap_ra = (map_ra - map_ra[map_ra.size / 2]) * np.pi / 180.
    map_ra = (map_ra - map_ra[map_ra.size / 2]) * np.pi / 180.
    map_dec = (map_dec - map_dec[np.sqrt(map_dec.size) / 2]) * np.pi / 180.
    # Interpolate map onto simulation cube grid
    gridded_map = interp.griddata(np.array([map_ra, map_dec]).T, map_data,
                                  np.array([sim_ra, sim_dec]).T,
                                  method=interp_method)
    # Scale the interpolation outputs to simulation frequencies
    interp_gridded_map = freq_interp(gridded_map, sim_freq)
    # Write an output file
    out = np.array([sim_ra, sim_dec, sim_freq, sim_data,
                   interp_gridded_map, sim_data + interp_gridded_map]).T
    np.savetxt(outfile, out,
               fmt='%15.12f %15.12f %12.7f %15.12f %17.12f %17.12f')


def freq_interp(base_data, interp_freq, base_freq=408., index=2.5):
    """Interpolate base_data at base_freq to interp_freq
    using power law of -index
    """
    return base_data * (interp_freq / base_freq) ** (-index)