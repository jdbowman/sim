"""
Store system and configuration information for running MAPS.

"""

import os

# Telescopes info
MWA128_LAT = -26.7012706
MWA128_LON = 116.670456
MWA128_EL = 365.76
VLA_LAT = 252.3210278
VLA_LON = 34.025778
VLA_EL = 2125.3704

# MAPS config
HOME = os.environ['HOME']
MAPS_DIR = HOME + '/local/src/MAPS'
ARRAY_DIR = MAPS_DIR + '/array'
ARRAY_CONFIG = {'MWA_128': ARRAY_DIR + '/mwa_128_crossdipole_gp_20110225.txt',
                'VLA_D': ARRAY_DIR + '/VLA_D.txt'}
ARRAY_LOC = {'MWA_128': (str(MWA128_LAT), str(MWA128_LON), str(MWA128_EL)),
             'VLA_D': (str(VLA_LAT), str(VLA_LON), str(VLA_EL))}
GHA = {'MWA_128': -7.778038703703704,
       'VLA_D': 16.821401851851853}
