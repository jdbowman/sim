"""
This module provides and (probably) complete and automate object to work on
drift scan simulations of the MWA using MAPS. It requires a python wrapper for MAPS (maps.py).

"""
import os
import astro
import settings as s
import numpy as np
import maps
import multiprocessing
from datetime import datetime


class _InputError(Exception):
    """Class for input error exceptions.

    """
    def __init__(self, message, func_name, object_name):
        self.err = message
        self.func_name = func_name
        self.object_name = object_name

    def __str__(self):
        return '{0}: {1}.{2}:\n>>>> {3}'\
               .format(str(datetime.now()), self.object_name, self.func_name,
                       self.err)


class Drift:
    """
    This class store information of a MAPS drift scan simulation of the MWA 128T
    (functionality can be expanded later) and a wrapper of maps.py to be use
    as a "single click" MAPS simulation

    """
    def __init__(self, target_ra, target_ha, sky_img=None, oobs=None,
                 pointing_center='zenith', fov_size=(412530.0, 412530.0),
                 duration=2.0, frequency=140.0, corr_int_time=1.0,
                 corr_chan_bw=0.04, scan_start='gha', site='MWA_128', name=None):
        """
        Initialize a drift scan

        Parameters
        ----------
        sky_img: string
            Name and path of the input sky image. The image must be a SIN
            projection of equal size in FITS format
        target_ra: float
            Right ascension at the center of a star/target to be observed
            [decimal hours]
        target_ha: float
            Hour angle of the star/target from the zenith meridian
            at the star of an observation [decimal hours]
        pointing_center: 'zenith' or ('hh:mm:ss', 'dd:mm:ss') (NOT IMPLEMENTED)
            Pointing center of the drift scan. Default mode 'zenith' uses
            the zenith calculated from target_ra and target_ha as a pointing.
        fov_size: array-like of float, 2 elements, optional
            Angular size of the field of view (RA[arcseconds, Dec[arcseconds])
            This is equivalent to the angular size of the sky image.
            The default value is exactly 2 radians thus equal to whole sky
            in SIN projected image.
        duration: float
            Duration of an observation [second]
        frequency: float
            Observing frequency [MHz]
        corr_int_time: float
            Correlator integration time [second]
        corr_chan_bw: float
            Correlator channel bandwidth [kHz]
        scan_start: 'gha' or 'year:day-of-year:hour:minute:second'
            Start time of a scan. If 'gha', the program will determine a
            Greenwich hour angle from the site keyword and use GHA convention
            to use relative time hard-coded in visgen.
            Else, user can give an absolute by giving a string of
            'year:day-of-year:hour:minute:second'. When specifying the absolute
            time, visgen will calculate the HA from this using the array
            location, but precision effects such as precession and nutation
            are not included.
        site: 'MWA_128', 'VLA_D'
            Defining array to use.
        name: string, optional
            Name of the observation. The name of sky_img - ".fits" with out path
            will be use if None

        """
        self.site = site
        if pointing_center == 'zenith':
            self.FOV_center_RA = astro.h2hms(target_ra + target_ha)
            self.FOV_center_Dec = astro.d2dms(float(s.ARRAY_LOC[self.site][0]))
        else:
            self.FOV_center_RA = pointing_center[0]
            self.FOV_center_Dec = pointing_center[1]
        self.FOV_size_RA = str(fov_size[0])
        self.FOV_size_Dec = str(fov_size[1])
        self.Frequency = str(frequency)
        self.Corr_int_time = str(corr_int_time)
        self.Corr_chan_bw = str(corr_chan_bw)
        self.Channel = self.Frequency + ':' + self.Corr_chan_bw
        if scan_start == 'gha':
            self.Scan_start = 'GHA {0:f}'.format(s.GHA[site])
        else:
            self.Scan_start = scan_start
        self.Scan_duration = str(duration)
        self.Time_cells = '0'
        self.Freq_cells = '0'
        if name is None:
            self.name = sky_img.rsplit('/', 1)[-1][0:-5]
        else:
            self.name = name
        self.sky_img = sky_img
        self.oobs = oobs
        self.spec_file = None
        self.vis_in = None
        self.vis_out = None
        self.vislog = None
        self.uvfits = None

        self.__spec = ''
        self.update_spec()
        self.__log = ''

    def __str__(self):
        return self.__spec + self.__log

    def update_spec(self):
        speckeys = np.array(['FOV_center_RA', 'FOV_center_Dec', 'FOV_size_RA',
                             'FOV_size_Dec', 'Corr_int_time', 'Corr_chan_bw',
                             'Scan_start', 'Scan_duration', 'Channel'])
        header = ('# {0}\n'
                  '# MAPS drift scan simulation\n'
                  '# name: {1}\n'
                  '# sky image: {2}\n'
                  '# OOB sources: {3}\n'
                  '# sky uvgrid: {4}\n'
                  '# visibility: {5}\n'
                  '# uvfits: {6}\n'
                  '# visgen log: {7}\n'
                  '# visgen specification file: {8}\n'
                  .format(str(datetime.now()), self.name, self.sky_img,
                          self.oobs, self.vis_in, self.vis_out, self.uvfits,
                          self.vislog, self.spec_file))
        spec = ''
        for k in speckeys:
            spec += '{0} = {1}\n'.format(k, self.__dict__.get(k))
        self.__spec = header + spec + 'Endscan\n\n'

    def print_spec(self):
        print self.__spec

    def write_spec(self):
        """
        Write observation specification file (*.ospec)

        Parameters
        ----------
        filename: str, optional
            Name of the ospec file
        verbose: boolean, optional
            Print the ospec file if True

        """
        outfile = self.name + '.ospec'
        self.spec_file = outfile
        self.update_spec()
        self.append_log('# $> write_spec()\n'
                        '# >>> visgen spec: {0}\n'
                        .format(self.spec_file))
        with open(outfile, 'w') as f:
            f.write(self.__spec)

    def append_log(self, string):
        self.__log += '# {0}\n{1}\n'.format(str(datetime.now()), string)

    def print_log(self):
        print self.__str__()

    def write_log(self):
        self.update_spec()
        self.append_log('# $> write_log()\n'
                        '# >>> log file: {0}\n'
                        .format(self.name + '.log'))
        with open(self.name + '.log', 'w') as f:
            f.write(self.__str__())

    def im2uv(self):
        if self.sky_img is None:
            raise _InputError('No imagae file', self.im2uv.__name__,
                              self.__name__)
        else:
            print '# im2uv: ' + self.name
            maps.im2uv(self.sky_img, verbose=False)
            self.vis_in = self.sky_img.rsplit('/', 1)[-1][0:-5] + '.dat'
            self.update_spec()
            self.append_log('# $> im2uv({0})\n'
                            '# >>> sky uvgrid: {1}\n'
                            .format(self.sky_img, self.vis_in))

    def visgen(self, mpi=1):
        if self.spec_file is None:
            raise _InputError('No oobs file')
        if self.vis_in is None and self.oobs is None:
            raise _InputError('Neither uvgrid file nor oob source list exist.',
                              self.visgen.__name__, self.name)
        else:
            print '# visgen: ' + self.name
            maps.visgen(self.name, self.spec_file, oobs=self.oobs,
                        uvgrid=self.vis_in, mpi=mpi, site=self.site)
            self.vis_out = self.name + '.vis'
            self.vislog = self.name + '.vislog'
            self.update_spec()
            self.append_log('# $> visgen()\n'
                            '# >>>> visgen visibility: {0}\n'
                            '# >>>> visgen log file: {1}\n'
                            .format(self.vis_out, self.vislog))

    def maps2uvfits(self):
        if self.vis_out is None:
            raise _InputError('visibility from visgen is not present',
                              self.maps2uvfits.__name__, self.name)
        else:
            print '# maps2uvfits: ' + self.name
            maps.maps2uvfits(self.vis_out, site=self.site, verbose=False)
            self.uvfits = self.name + '.uvfits'
            self.update_spec()
            self.append_log('# $> maps2uvfits({0})\n'
                            '# >>>> uvfits: {1}\n'
                            .format(self.vis_out, self.uvfits))

    def go(self):
        if self.sky_img is not None:
            self.im2uv()
        self.write_spec()
        self.visgen()
        if self.vis_in is not None:
            os.remove(self.vis_in)
            self.append_log('# remove ' + self.vis_in)
        self.maps2uvfits()
        self.write_spec()
        self.write_log()


def __call_go(instance):
    return instance.go()


def batch_drift(instance, nthreads=4):
    pool = multiprocessing.Pool(nthreads)
    pool.map(__call_go, instance)
    pool.close()
    pool.join()
