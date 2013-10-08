import os
import numpy as np
from hpm2sin import hpm2sin
freqs = [142.735 + i * .04 for i in range(36)]
hpmfile = ['/data3/piyanat/model/21cm/hpm/delta21cm_ComovingGrid_11-Aug-2013_hpx_{:.3f}.fits'.format(f) for f in freqs]

fov_ra = 4.0
fov_dec = -26.7012706
ha = np.arange(-40 * 94, 40 * 94, 94) / (60. * 60)                
ra = fov_ra + ha
dec = np.ones_like(ra) * fov_dec
multiplier = np.ones(dec.size)
for hpm, f in zip(hpmfile, freqs):
    outdir = '/data3/piyanat/model/21cm/sin/{:.3f}/'.format(f)
    if not os.path.exists(outdir):
        os.makedirs(outdir)    
    fitsfile = ['{:s}delta21cm_ComovingGrid_11-Aug-2013_{:.3f}_{:.3f}.fits'.format(outdir, h, f) for h in ha]
    hpm2sin(hpm, fitsfile, ra, dec, nthreads=8)
