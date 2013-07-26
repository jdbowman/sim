"""
A script template to run multiple drift scan.

"""
import os
import drift
import numpy as np
from glob import glob

frequencies = [i * 1.28 - 0.64 + (4 - 1) * 0.005 for i in range(63, 157)]
ra = 4.0
hour_angles = np.arange(-40 * 94, 40 * 94, 94) / (60. * 60)
WORK_DIR = '/home/piyanat/workspace/run/zenith_drift/'
IMG_DIR = '/data2/piyanat/diffuse/sin/haslam/'

for freq in frequencies:
    print '# ==================================================================='
    print '# Simulting at {:.3f} MHz'.format(freq)
    print '# -------------------------------------------------------------------'
    images = glob(IMG_DIR + '{:.3f}/*.fits'.format(freq))
    images.sort()
    names = [img.rsplit('/', 1)[-1][:-5] for img in images]
    OUT_DIR = WORK_DIR + '{:.3f}/'.format(freq)
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
        os.chdir(OUT_DIR)
    print '# ' + os.getcwd()
    snapshot = [drift.Drift(ra, ha, sky_img=img, oobs=None, frequency=freq,
                            corr_int_time=2.0, corr_chan_bw=0.04, name=n)
                for ha, img, n in zip(hour_angles, images, names)]
    drift.batch_drift(snapshot, nthreads=8)
    os.chdir(WORK_DIR)
    print '# ' + os.getcwd()
