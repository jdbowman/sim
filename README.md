MAPS and tools for redshifted 21 cm instrument simulation

DIRECTORIES
-----------
MAPS: MIT Array Performace Simulator. Use to simulate observation of the MWA
    and other telescope. The most recent version of MAPS can be check out from
    the Haystack server

    svn checkout svn+ssh://mwa-lfd.haystack.mit.edu/svn/MAPS

simsky: consists of srcgen, galmask and galfix. Can be used to simulate
    foreground to feed into MAPS

pymodule: python modules, consisting of python wrapper of MAPS tasks and tools
    to do other stuff

scripts: template of python scripts that utilize pymodule and scripts to do
    other stuff such as extrapolating the Haslam map. Some scripts are
    executable expected in pymodule, so you should add this directory to your
    PATH
