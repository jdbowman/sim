import numpy as np
import pyfits as pf
import pywcs as wcs
import scipy.interpolate as interp


def ReadSim(simfile, cube=False):
    """Read sky coordinate (RA, Dec), frequencies and temperature from 21 cm
    simulation cube.

    Parameters
    ----------
    simfile: str
        A simulation cube file
    cube: boolean, optional
        If true, outputs will be in cube (x, y, freq) formats.

    Returns
    -------
    out: ndarray
        RA, Dec, frequencies and temperatures read from the simulation cube
    """
    ra, dec, freq, temp = np.genfromtxt(simfile, unpack=True)
    if cube:
        return ra.reshape((128, 128, 128)), dec.reshape((128, 128, 128)), \
        freq.reshape((128, 128, 128)), temp.reshape((128, 128, 128))
    else:
        return ra, dec, freq, temp


def ReadMap(mapfile):
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


def AddMap2Sim(mapfile, simfile, outfile, interp_method='nearest'):
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
    map_ra, map_dec, map_data = ReadMap(mapfile)
    sim_ra, sim_dec, sim_freq, sim_data = ReadSim(simfile)
    # Convert RA and Dec of the diffuse map into a cube-centric systemmap_ra = (map_ra - map_ra[map_ra.size / 2]) * np.pi / 180.
    map_ra = (map_ra - map_ra[map_ra.size / 2]) * np.pi / 180.
    map_dec = (map_dec - map_dec[np.sqrt(map_dec.size) / 2]) * np.pi / 180.
    # Interpolate map onto simulation cube grid
    gridded_map = interp.griddata(np.array([map_ra, map_dec]).T, map_data,
                                  np.array([sim_ra, sim_dec]).T,
                                  method=interp_method)
    # Scale the interpolation outputs to simulation frequencies
    interp_gridded_map = FreqInterp(gridded_map, sim_freq)
    # Write an output file
    out = np.array([sim_ra, sim_dec, sim_freq, sim_data,
                   interp_gridded_map, sim_data + interp_gridded_map]).T
    np.savetxt(outfile, out,
               fmt='%15.12f %15.12f %12.7f %15.12f %17.12f %17.12f')


def FreqInterp(base_data, interp_freq, base_freq=408., index=2.5):
    """Interpolate base_data at base_freq to interp_freq
    using power law of -index
    """
    return base_data * (interp_freq / base_freq) ** (-index)


# def Smooth(sky, xpoints=128, ypoints=128):
#   x = np.arange(sky.shape[0])
#   y = np.arange(sky.shape[1])
#   z = sky.reshape(sky.size)
#   f = interp.interp2d(x, y, z, kind='cubic')
#   xnew = np.linspace(np.min(x), np.max(x), xpoints)
#   ynew = np.linspace(np.min(y), np.max(y), ypoints)
#   return f(xnew, ynew)


# def Run():
#   basesky = ReadMap()
#   ra, dec, freq, temp = ReadSim('delta_21cm_1gpch_z6.90_ng128_AAFT.txt')
#   newsky = np.zeros((basesky.shape[0], basesky.shape[1], freq.size))
#   smoothsky = np.zeros((temp.shape[0], temp.shape[1], freq.size))
#   for i in range(len(freq)):
#       newsky[:, :, i] = FreqInterp(basesky, freq[i])
#       smoothsky[:, :, i] = Smooth(newsky[:, :, i])
#   return smoothsky


# def Grid1():
#   ra, dec, freqgrid, temp = ReadSim('Corrected21cmSignal_AAFT.txt')
#   freq = freqgrid[0, 0, :]
#   sky = ReadMap()
#   interp_sky = np.zeros((sky.shape[0], sky.shape[1], 128))
#   newsky = np.zeros((128, 128, 128))
#   for i in range(128):
#       interp_sky[:, :, i] = FreqInterp(sky, freq[i])
#   x, y = np.mgrid[0:sky.shape[0], 0.:sky.shape[1]]
#   x = ((x - (sky.shape[0] / 2. - 1)) * 10. / sky.shape[0]) * (np.pi / 180.)
#   y = ((y - (sky.shape[1] / 2. - 1.)) * 10. / sky.shape[1]) * (np.pi / 180.)
#   points = np.transpose(np.array([x.reshape(sky.size), y.reshape(sky.size)]))
#   gridx = np.zeros((128 * 128, 128))
#   gridy = np.zeros((128 * 128, 128))
#   for i in range(128):
#       gridx[:, i] = ra[:, :, i].reshape(128 * 128)
#       gridy[:, i] = dec[:, :, i].reshape(128 * 128)
#       newsky[:, :, i] = interp.griddata(points,
#                       interp_sky[:, :, i].reshape(sky.size),
#                       (gridx[:, i], gridy[:, i]),
#                       method='linear').reshape(128, 128)
#   return ra, dec, freqgrid, temp, newsky


# def Grid2():
#   ra, dec, freqgrid, temp = ReadSim('Corrected21cmSignal_AAFT.txt')
#   freq = freqgrid[0, 0, :]
#   sky = ReadMap()
#   interp_sky = FreqInterp(sky, freq[0])
#   newsky = np.zeros((128, 128, 128))
#   x, y = np.mgrid[0:sky.shape[0], 0.:sky.shape[1]]
#   x = ((x - (sky.shape[0] / 2. - 1)) * 10. / sky.shape[0]) * (np.pi / 180.)
#   y = ((y - (sky.shape[1] / 2. - 1.)) * 10. / sky.shape[1]) * (np.pi / 180.)
#   points = np.transpose(np.array([x.reshape(sky.size), y.reshape(sky.size)]))
#   gridx = ra[:, :, 0].reshape(128 * 128)
#   gridy = dec[:, :, 0].reshape(128 * 128)
#   newsky[:, :, 0] = interp.griddata(points, interp_sky.reshape(sky.size),
#                   (gridx, gridy), method='linear').reshape(128, 128)
#   for i in range(1, 128):
#       newsky[:, :, i] = FreqInterp(newsky[:, :, 0], freq[i], freq[0])
#   return ra, dec, freqgrid, temp, newsky


# def Combine(name, ra, dec, freq, temp, sky):
#   out = np.transpose(np.array([ra, dec, freq, temp, sky, temp + sky])\
#       .reshape(6, 128 * 128 * 128))
#   np.savetxt(name, out, fmt='%0.6f\t%0.6f\t%3.2f\t%0.8f\t%3.8f\t%3.8f')
