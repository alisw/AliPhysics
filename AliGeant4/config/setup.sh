# $Id$
# ----------------------------------------------------------------
# This script sets the default environment variables for
# Alice Geant4 based prototype
#   
# by I. Hrivnacova, 18.8.1998
#
# sh version modified by I. Gonzalez 2.3.2000
#
# modified by fca Sep 2002

export MCINSTALL=$ALICE/geant4vmc

# This is architecture dependent...
SYSTEM=`uname`

#
# path to shared libraries
# 
if [ "$SYSTEM" = "HP-UX" ]; then
  SHLIBVAR=$SHLIB_PATH
  SHLIBVARNAME=SHLIB_PATH
fi  
if [ "$SYSTEM" = "Linux" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi
if [ "$SYSTEM" = "OSF1" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi
if [ "$SYSTEM" = "SunOS" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi

if [ "`echo ${SHLIBVAR} | grep ${G4INSTALL}/lib/${G4SYSTEM} `" = "" ]; then
  SHLIBVAR="${SHLIBVAR}:${G4INSTALL}/lib/${G4SYSTEM}"
fi
if [ "`echo ${SHLIBVAR} | grep ${MCINSTALL}/lib/tgt_${G4SYSTEM} `" = "" ]; then
  SHLIBVAR="${SHLIBVAR}:${MCINSTALL}/lib/tgt_${G4SYSTEM}"
fi
if [ "`echo ${SHLIBVAR} | grep ${CLHEP_BASE_DIR}/lib `" = "" ]; then
  SHLIBVAR="${SHLIBVAR}:${CLHEP_BASE_DIR}/lib"
fi

export $SHLIBVARNAME=$SHLIBVAR

#
# Remove unneeded variables. If this is not done the vars. remain in the env.
#

unset SYSTEM
unset SHLIBVAR
unset SHLIBVARNAME
