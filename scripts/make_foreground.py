"""Make diffuse foreground by extrapolating the Haslam map.
This script is now packed into a set of tools to run MAPS.
"""
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

hpm_in_name = 'lambda_haslam408_dsds.fits'
hpm_out_name = 'haslam_xpt_N4096.fits'
nside_out = 4096                                                                                               
hpm_in_lmax = None
fitrange = (36, 90)

# Read the input map
hpm = hp.read_map(hpm_in_name)

# Calculate the map power spectrum
hpm_cl, hpm_alm = hp.anafast(hpm, lmax=hpm_in_lmax, alm=True)

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
hpm_extrap_cl = hp.anafast(hpm_extrap)
plt.loglog(np.arange(hpm_extrap_cl.size), hpm_extrap_cl, 'y', label='Model')
plt.loglog(np.arange(hpm_cl.size), hpm_cl, 'b--', lable='Haslam')
plt.loglog(np.arange(diff_cl.size), diff_cl, 'r-.', lable='Small Structure')
plt.legend(loc='upper right')
plt.axvline(fitrange(0), c='k', ls=':')
plt.axvline(fitrange(1), c='k', ls=':')
hp.mollview(hpm_original, title='Haslam Original N4096')
hp.mollview(hpm_diff, title='Small Structure')
hp.mollview(hpm_extrap, title='Model = Haslam + Small Struture')
