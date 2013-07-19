"""scale extracted haslam maps

"""
import os
import numpy as np
from img import freq_interp_fits
from glob import glob

SIN_DIR = '/home/piyanat/workspace/diffuse/sin/'
GAMMA_DIR = SIN_DIR + 'gamma/'
BETA_DIR = SIN_DIR + 'beta/'
HASLAM_DIR = SIN_DIR + 'haslam/408/'

fov_ra = 4.0
fov_dec = -26.7012706
ha = np.arange(-40 * 94, 40 * 94, 94) / (60. * 60)
ra = fov_ra + ha
dec = np.ones_like(ra) * fov_dec

# freqs = [i * 1.28 - 0.64 + (4 - 1) * 0.005 for i in range(63, 157)]
# fine channels over a continuous 8 MHz bandwidth. We already interpolate to
# the first (126.095 MHz) and the last (133.775 MHz), so we skip those.
freqs = [126.095 + (i * 0.04) for i in range(193)][1:-2]

# get file list
haslam_sin = np.array(glob(HASLAM_DIR + '*'))
haslam_sin.sort()
beta_sin = np.array(glob(BETA_DIR + '*'))
beta_sin.sort()
gamma_sin = np.array(glob(GAMMA_DIR + '*'))
gamma_sin.sort()

# Check existing directory
for f in freqs:
    outdir = '{0}haslam/{1:.3f}/'.format(SIN_DIR, f)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

for h, b, g in zip(haslam_sin, beta_sin, gamma_sin):
    outfits = ['{0}haslam/{1:.3f}/{2}_{1:.3f}.fits'
               .format(SIN_DIR, f, h.rsplit('/', 1)[-1][:-5]) for f in freqs]
    freq_interp_fits(h, 408.0, freqs, outfits=outfits, beta=b,
                     gamma=g, T2Jy=True)