"""
Python wrapper of main tasks on the MIT Array Performance Simylator (MAPS)

"""
import settings as s
from subprocess import Popen, PIPE, call, STDOUT


def _save_string(filename, string):
    with open(filename, 'w') as f:
        f.write(string)


class _Error(Exception):
    """Base class for error exceptions in this module.

    """
    pass


class _InputError(_Error):
    """Class for input error exceptions.

    """
    pass


class _VisgenError(_Error):
    """Class for handling error from visgen

    """
    def __init__(self, message, filename):
        self.err = message
        self.errfile = filename

    def __str__(self):
        return 'The following message is saved to {0} \n {1}'\
            .format(self.errfile, self.err)


def im2uv(fitsfile, vis=None, multiplier=None, padzeropixels=None,
          verbose=True):
    """Convert a FITS image into a visibility grid format appropriated for
    visgen input via MAPS_im2uv.

    Parameters
    ----------
    fitsfile: str
        Name of an input FITS image.
    vis: str, optional
        Name of an output visibility grid file (*.dat).
        If None, vis = (fitsfile - '.fits') + '.dat'
    multiplier: float, optional
        Multiply the pixel value of the input FITS image by this value.
        Can be used to adjust to unit to Jy/sr as expected by visgen
    padzeropixels: int, optional
        Pad the input FITS image with a given number of zero pixels on each side
    verbose: boolean, optional
        Set verbosity of the program.
        If True, stout and stderr will be dumped onto terminal in real time.
        No log file will be save.
        If False, no terminal dump. All stdout and stderr is save to a file
        named (vis - '.dat') + .im2uvlog

    """
    if vis is None:
        vis = fitsfile.rsplit('/', 1)[-1][0:-5] + '.dat'
    cmd = ['MAPS_im2uv', '-i', fitsfile, '-o', vis]
    if multiplier is not None:
        cmd += ['-n', str(multiplier)]
    if padzeropixels is not None:
        cmd += ['-p', str(padzeropixels)]
    if verbose:
        call(cmd)
    else:
        run = Popen(cmd, stdout=PIPE, stderr=STDOUT)
        stdout = run.communicate()[0]
        if stdout != '':
            logfile = vis.rsplit('/', 1)[-1][0:-4] + '.im2uvlog'
            _save_string(logfile, stdout)


def maps2uvfits(vis, uvfits=None, site='MWA_128', arrayloc=None, arrayconf=None,
                verbose=True):
    """Convert visgen visibility grid to AIPS uvfits via maps2uvfits

    """
    if arrayconf is None:
        arrayconf = s.ARRAY_CONFIG[site]
    if arrayloc is None:
        arrayloc = s.ARRAY_LOC[site]
    if uvfits is None:
        uvfits = vis.rsplit('/', 1)[-1][0:-4] + '.uvfits'
    cmd = ['maps2uvfits', vis, uvfits, arrayloc[0], arrayloc[1], arrayloc[2],
           arrayconf]
    if verbose:
        call(cmd)
    else:
        run = Popen(cmd, stdout=PIPE, stderr=STDOUT)
        stdout, stderr = run.communicate()
        if stdout != '':
            logfile = vis.rsplit('/', 1)[-1][0:-4] + '.maps2uvfitslog'
            _save_string(logfile, stdout)


def visgen_base(args):
    """This function is a base for visgen. It uses subprocess to call visgen
    from command-line passing cmd argument

    """
    prefix, cmd = args
    run = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = run.communicate()
    _save_string(prefix + '.vislog', stdout)
    if stderr != '':
        _save_string(prefix + '.viserr', stderr)
        raise _VisgenError(stderr, prefix + '.viserr')


def visgen(prefix, spec, oobs=None, uvgrid=None, mpi=1, site='MWA_128'):
    """Wrapper of visgen_base. This allow running multiple visgen jobs in
    parallel on a single node via multiprocessing module., i.e. utilizing
    multi-core processing instead of a for loop.

    prefix: string
        prefix of the output file
    spec: string
        name and path of visgen specifivation file (ospec file)
    oobs: string
        name and path of out-of-bound source list
    uvgrid: string
        name and path of the input uvgrid (dat file) produced by im2uv or brigen
    mpi: integer, optional
        if > 1 visgen will be run with mpi

    """
    arrayconf = s.ARRAY_CONFIG[site]
    case = {'oobs_only': oobs is not None and uvgrid is None,
            'uvgrid_only': oobs is None and uvgrid is not None,
            'both': oobs is not None and uvgrid is not None}
    if case['oobs_only']:
        cmd = ['visgen', '-n', prefix, '-s', site, '-A', arrayconf,
               '-V', spec, '-O', oobs, '-Z', '-N', '-m', '0']
    elif case['uvgrid_only']:
        cmd = ['visgen', '-n', prefix, '-s', site, '-A', arrayconf,
               '-V', spec, '-G', uvgrid, '-N', '-m', '0']
    elif case['both']:
        cmd = ['visgen', '-n', prefix, '-s', site, '-A', arrayconf,
               '-V', spec, '-O', oobs, '-G', uvgrid, '-N', '-m', '0']
    else:
        raise _InputError('Need at least point sources file or uvgrid \
                          foreground input')
    if mpi > 1:
        cmd = ['mpirun', '-n', str(mpi)] + cmd
    args = [prefix, cmd]
    visgen_base(args)
