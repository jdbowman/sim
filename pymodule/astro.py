"""Astronomical tools

"""

import datetime


def h2hms(h):
    """
    Convert decimal hours to hh:mm:ss

    """
    inverse = ''
    if h < 0:
        inverse = '-'
    hms = str(datetime.timedelta(hours=abs(h))).rsplit(', ')[-1]
    return inverse + hms


def d2dms(d, delimiter=':', precision=6):
    """
    Convert decimal degrees to dd:mm:ss

    """
    inverse = ''
    if d < 0:
        inverse = '-'
    minutes, seconds = divmod(abs(d) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    return '{0:s}{1:.0f}{4}{2:02.0f}{4}{3:0{5:d}.{6:d}f}'\
        .format(inverse, degrees, minutes, seconds, delimiter, 3 + precision,
                precision)


def lst2gha(lst, site_long=116.670456):
    """
    Convert LST in decimal degree to GHA in decimal hours. Longitude of
    the observer can be specify with site_long = +/- decimal degree, where
    + is east and - is west of the prime meridian. Else assume MWA 128T
    location, site_long = 116.670456.

    """
    gha = lst/15. - site_long/15.
    if gha > 24.0:
        gha -= 24.0
    elif gha < 0.0:
        gha += 24.0
    return gha
