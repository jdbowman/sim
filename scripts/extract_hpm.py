import os
import numpy as np
from img import freq_interp_fits, extract_hpm
from glob import glob

SOURCE_DIR = 'fullsky/'
SIN_DIR = 'sin/'
GAMMA_DIR = SIN_DIR + 'gamma/'
BETA_DIR = SIN_DIR + 'beta/'
HASLAM_DIR = SIN_DIR + 'haslam/408/'
haslam = 'haslam_xtp_N4096.fits'
beta = 'beta_0.04sigma_N4096.fits'
gamma = 'gamma_0.01sigma_N4096.fits'

fov_ra = 4.0
fov_dec = -26.7012706
ha = np.arange(-40 * 94, 40 * 94, 94) / (60. * 60)
ra = fov_ra + ha
dec = np.ones_like(ra) * fov_dec

# freqs = [i * 1.28 - 0.64 + (4 - 1) * 0.005 for i in range(63, 157)]

# extracting
for files, dirs in zip([haslam, beta, gamma], [HASLAM_DIR, BETA_DIR, GAMMA_DIR]):
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    outfiles = ['{0}{1}_sin_zd_eor1_{2:.3f}.fits'
                .format(dirs, files.rsplit('_', 1)[0], h) for h in ha]
    extract_hpm(SOURCE_DIR + files, outfiles, ra, dec, nthreads=8)

# scale extracted haslam maps
# haslam_sin = np.array(glob(HASLAM_DIR + '*'))
# haslam_sin.sort()
# beta_sin = np.array(glob(BETA_DIR + '*'))
# beta_sin.sort()
# gamma_sin = np.array(glob(GAMMA_DIR + '*'))
# gamma_sin.sort()
# for f in freqs:
#     path = '{0}haslam/{1:.3f}/'.format(SIN_DIR, f)
#     if not os.path.exists(path):
#         os.makedirs(path)
#     for h, b, g in zip(haslam_sin, beta_sin, gamma_sin):
#         p = path + h.rsplit('.', 1)[0]
#         print p
