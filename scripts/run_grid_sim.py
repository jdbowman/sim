from multiprocessing import Pool
from grid_sim import grid_sim

freqs = [140.175 + i * .04 for i in range(200)]
fitsfile = ['/home/piyanat/workspace/model/21cm/healpix/delta21cm_ComovingGrid_11-Aug-2013_hpx_{:.3f}.fits'.format(f) for f in freqs]
simfile = ['/home/piyanat/workspace/model/21cm/cube/delta21cm_ComovingGrid_11-Aug-2013.npy' for i in range(len(freqs))]

def call_grid_sim(args):
    print args
    simfile, fitsfile, freq = args
    grid_sim(simfile, fitsfile, freq)

workers = Pool(12)
workers.map(call_grid_sim, zip(simfile, fitsfile, freqs))
