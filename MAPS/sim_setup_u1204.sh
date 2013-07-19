#
# sim_setup.sh
#
# MAPS setup script
#
# Created: Stone Age (petroglyphs)
#
# Modified:
# 13 Jan 2011 by L. Benkevitch. Added optional cmdline parameter, the
#     MAPS directory with working copy mame. Users are advised to 
#     include in their '.bashrc' file a line
#
#     source sim_setup.sh MAPS
#
#     where MAPS is either MAPS or their MAPS working copy name (e.g., maps,
#     maps_new, Maps, MAPS_test etc. I, for example, usually use just maps) 
#     If the MAPS parameter is absent, default is MAPS.
#     Also, commented out printing
#     "HAVE_CBLAS is set to ..."
#     "HAVE_WCSLIB is set to  ..."
#
#     P.S. In case of multiple MAPS working copies do not forget to 
#     source sim_setup.sh with the current MAPS name prior to compilation!
#     
# 19 Jul 2013 by Piyanat Kittiwisit. Create sim_setup_u1204.sh specifically
#     to be used in a system with Ubuntu 12.04 LTS. Environment variables are
#     set to utilize libraries installed using apt-get from the Unbuntu
#     repository. Users still need to edit the location of the MAPS directory.
#      

export SIM_CFLAGC="-O"
export SIM_CFLAGL=-O

if [ "$1" = "" ] ; then
    export MAPS_WCNAME='MAPS'
else
    export MAPS_WCNAME=$1
fi

# you MUST change the line below to point to the base
# directory of the MAPS installation
export SIM=/usr/global/sim/$MAPS_WCNAME
echo MAPS directory is set to '$SIM='\'$SIM\'

export HAVE_CBLAS="FALSE"
#export HAVE_WCSLIB="FALSE"
export HAVE_WCSLIB="TRUE"

# you MUST change the two lines below to point to the location of
# header files and libraries for wcslib
if [ "$HAVE_WCSLIB" = "TRUE" ] ; then
#  echo "HAVE_WCSLIB is set to TRUE"
  export WCS_INC=/usr/include/wcslib
  export WCS_LIB=/usr/lib
#else
#  echo "HAVE_WCSLIB is set to FALSE"
fi

# It seems that you do not need CBLAS to run MAPS
if [ "$HAVE_CBLAS" = "TRUE" ] ; then
#  echo "HAVE_CBLAS is set to TRUE"
  export CBLAS_LIB=
  export CBLAS_INC=
#else
#  echo "HAVE_CBLAS is set to FALSE"
fi

# This is the default location of cfitsio if you install it through apt-get.
# You will need to modify it if you install it other ways
export CFITSIO_INC=/usr/include
export CFITSIO_LIB=/usr/lib/x86_64-linux-gnu

# Other environment variables to make running MAPS easier
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
#if [ -d "$SIM_BIN" ] ; then
#    export PATH="$PATH:$SIM_BIN"
#fi

# Avoid path bloat
# Remove paths with MAPS/bin or maps/bin
PATH=$(echo $PATH | sed -e 's;:\?/home/benkev/maps/bin;;gI' -e 's;/home/benkev/maps/bin:\?;;gI')
#if [ ! `echo $PATH | grep $SIM_BIN` ] ; then
    export PATH="$PATH:$SIM_BIN"
#fi
