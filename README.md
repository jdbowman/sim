MAPS and tools for redshifted 21 cm instrument simulation

Directories
-----------
- **MAPS**:
MIT Array Performace Simulator. Use to simulate observation of the MWA
and other telescope. The most recent version of MAPS can be check out from
the Haystack server

        svn checkout svn+ssh://mwa-lfd.haystack.mit.edu/svn/MAPS

- **simsky**:
consists of srcgen, galmask and galfix. Can be used to simulate
foreground to feed into MAPS

- **pymodule**:
python modules, consisting of python wrapper of MAPS tasks and
tools to do other stuff

- **scripts**:
template of python scripts that utilize pymodule and scripts to do other stuff
such as extrapolating the Haslam map. Some scripts are executable expected in
pymodule, so you should add this directory to your PATH

Installation
------------
MAPS and simsky need to be complied separately. See README in each directory
for more detail. Note that MAPS contain a c++ adaptation of the Angelica's
Global Sky Model in MAPS/maps_makesky/gsm. This program requires database files
larger than 100MB which cannot be store on the GitHub server. You **MUST**
download these database files from https://www.dropbox.com/sh/ggsk91mikslwa2g/v3elY_tyox before compiling MAPS

pymodule and scrips can be directly added to you PYTHONPATH.