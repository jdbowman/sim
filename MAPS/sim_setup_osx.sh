#
# sim_setup_osx.sh
#
# MAPS setup script
#
# Creator: P Kittiwisit. M setup template for Mac OS X installation based
# on sim_setup.sh
#    

if [ "$1" = "" ] ; then
    export MAPS_WCNAME='MAPS'
else
    export MAPS_WCNAME=$1
fi

export SIM_CFLAGC="-O"
export SIM_CFLAGL=-O

# you MUST change the following part to point MAPS to location of
# dependent libraries
export WCS_LIB=/usr/local/lib
export WCS_INC=/usr/local/include/wcslib
export CBLAS_LIB=/usr/lib
export CBLAS_INC=/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
export CFITSIO_LIB=/usr/local/lib
export CFITSIO_INC=/usr/local/include
export NETCDF_LIB=/usr/local/lib
export NETCDF_INC=/usr/local/include


# you MUST change the line below to point to the base
# directory of the MAPS installation
# export SIM=/data/$USER/MAPS
export SIM=$HOME/local/src/$MAPS_WCNAME
echo MAPS directory is set to '$SIM='\'$SIM\'

export SIM_BIN=$SIM/bin
export SIM_LIB=$SIM/lib
export SIM_INC=$SIM/include
export VISDIR=$SIM/visdata

export TEXTDIR=$SIM/text
export ARRAYDIR=$SIM/array	
export LAYOUTDIR=$SIM/stn_layout


export GSMDIR=$SIM/maps_makesky/gsm
export FILE408=$SIM/maps_makesky/scripts/408POLN.txt


# Avoid path bloat
if [ -d "$SIM_BIN" ] ; then
    export PATH="$PATH:$SIM_BIN"
fi