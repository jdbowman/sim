#!/usr/bin/env python

"""
Program: grid_sim.py
Author: Piyanat Kittiwisit
        piyanat.kittiwisit@asu
Creatd: August 19, 2013

Required:
Tested on:

"""
import argparse
import numpy as np
import healpy as hp
from scipy.integrate import quad

def grid_sim(simfile, fitsfile, freq, nside=4096, simfile_type='npy',
             temp_field='temp', sep=',', col=3, sim_size=(128, 128, 128),
             sim_res=7.8125, coord_file=None):
    """
    Parameters
    ----------
    simfile: string
        Path and name of the file containing 21cm brightness temperature simulation cube
    fitsfile: string
        Name of the output HEALPix image
    freq: float
        frequency to the center of the shell of the sky of interest
    nside: integer
        NSIDE of the output HEALPix image. Must be a valid NSIDE for HEALPix.
    simfile_type: ascii or npy
        Type of the simulation file. Can be either an ascii table or a numpy binary file (npy).
        If ascii, one column must be brightness temperature sorted in z, y and then x
        npy file should contains numpy record array with one field as brightness temperature.
    temp_field: string
        Name of the brightness temperature field in the record array in the npy file
    sep: string
        Delimiter of the ascii table
    col: integer
        Column of the brightness temperature in the ascii file (zero base)
    sim_size: (integer, integer, integer)
        (x, y, z) sizes of the simulation cube
    sim_res: float
        Resolution of the simulation cube in Mpc/h'
    coord_file: string
       Save coordinate (xi, yi, zi) of the simulation cube mapped to HEALPix
       pixels to a numpy binary file

    """
    # Read simulation cube
    #temp = np.genfromtxt(args.sim_file, delimiter=args.sep, usecols=(args.col,))
    #temp.shape = args.sim_shape
    if simfile_type == 'npy':
        temp = np.load(simfile)[temp_field]
    elif simfile_type == 'ascii':
        temp = np.genfromtxt(simfile, delimiter=sep, usecols=(col,))
        temp.shape = sim_size

    # Determine the radial comoving distance r to the center of the shell
    # representing the frequency channel of interest
    H_0 = 70  # km / s / Mpc
    omega_m = 0.3
    omega_l = 1 - omega_m
    c = 299792458  # m / s
    f21 = 1420.40575177  # MHz
    z = f21 / freq - 1
    integrand = lambda z: 1 / np.sqrt(omega_m * (1 + z) ** 3 + omega_l)
    dc = (c / H_0) * quad(integrand, 0, z)[0] / 1000

    # Get the vector coordinates (x, y, z) of the HEALPIX pixels
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)
    x, y, z = hp.pix2vec(nside, pix)
    x *= dc
    y *= dc
    z *= dc

    # Determine the indexes (xi, yi, zi) of tiled simulation cube corresponding
    # to the vector coordinates (x, y, z) of the HEALPix pixels
    xi = np.mod(np.around(x / sim_res).astype(int), temp.shape[0])
    yi = np.mod(np.around(y / sim_res).astype(int), temp.shape[1])
    zi = np.mod(np.around(z / sim_res).astype(int), temp.shape[2])

    # Map the HEALPix pixels to the tiled cube
    out = temp[xi, yi, zi]
    hp.write_map(fitsfile, out, fits_IDL=False, dtype=np.float64)
    if coord_file is not None:
        np.save(coord_file, np.core.records.fromarrays([xi, yi, zi], names='xi,yi,zi'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    inparams = parser.add_argument_group('Input File Options')
    inparams.add_argument('simfile', type=str,
        help='Path and name of the file containing 21cm brightness temperature simulation cube')
    inparams.add_argument('--simfile_type', type=str, choices=['ascii', 'npy'], default='npy',
        help='Type of the simulation file. Can be either an ascii table or a numpy binary file (npy).'
             'If ascii, one column must be brightness temperature sorted in z, y and then x.'
             'npy file should contains numpy record array with one field as brightness temperature.')
    inparams.add_argument('--temp_field', type=str, metavar='field_name', default='temp',
        help='Name of the brightness temperature field in the record array in the npy file')
    inparams.add_argument('--sep', type=str, default=',', metavar='delimiter_string',
        help='Delimiter of the ascii table')
    inparams.add_argument('--col', type=int, default=3, metavar='column_number',
        help='Column of the brightness temperature in the ascii file (zero base)')
    inparams.add_argument('--sim_size', type=int, nargs=3, default='128 128 128'.split(),
        metavar=('sim_xsize', 'sim_ysize', 'sim_zsize'),
        help='Size of the simulation cube')
    inparams.add_argument('--sim_res', type=float, default=7.8125, metavar='resolution',
        help='Resolution of the simulation cube in Mpc/h')
    outparams = parser.add_argument_group('Output Parameters')
    outparams.add_argument('fitsfile', type=str,
       help='Name of the output HEALPix image')
    outparams.add_argument('freq', type=float,
       help='Frequency to the center of the shell of the sky of interest')
    outparams.add_argument('--nside', type=int, default=4096,
       help='NSIDE of the output HEALPix image')
    outparams.add_argument('--coord_file', type=str, metavar='filename',
       help='Save coordinate (xi, yi, zi) of the simulation cube mapped to HEALPix pixels')
    args = parser.parse_args()
    grid_sim(args.simfile, args.fitsfile, args.freq, nside=args.nside,
             simfile_type=args.simfile_type, temp_field=args.temp_field,
             sep=args.sep, col=args.col, sim_size=args.sim_size,
             sim_res=args.sim_res, coord_file=args.coord_file)
