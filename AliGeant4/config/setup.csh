# ----------------------------------------------------------------
# This script sets the default environment variables for
# Alice Geant4 based prototype. Based on the sh version.
# I. Gonzalez 03.09.2002

setenv MCINSTALL $ALICE/geant4vmc
setenv SYSTEM `uname`

#
# path to shared libraries
# 
if ( "$SYSTEM" == "HP-UX" ) then
  set SHLIBVAR=$SHLIB_PATH
  set SHLIBVARNAME=SHLIB_PATH
endif
if ( "$SYSTEM" == "Linux" ) then
  set SHLIBVAR=$LD_LIBRARY_PATH
  set SHLIBVARNAME=LD_LIBRARY_PATH
endif
if ( "$SYSTEM" == "OSF1" ) then
  set SHLIBVAR=$LD_LIBRARY_PATH
  set SHLIBVARNAME=LD_LIBRARY_PATH
endif
if ( "$SYSTEM" == "SunOS" ) then
  set SHLIBVAR=$LD_LIBRARY_PATH
  set SHLIBVARNAME=LD_LIBRARY_PATH
endif

if ( "`echo ${SHLIBVAR} | grep ${G4INSTALL}/lib/${G4SYSTEM} `" == "" ) then
  set SHLIBVAR="${SHLIBVAR}:${G4INSTALL}/lib/${G4SYSTEM}"
endif
if ( "`echo ${SHLIBVAR} | grep ${MCINSTALL}/lib/tgt_${G4SYSTEM} `" == "" ) then
  set SHLIBVAR="${SHLIBVAR}:${MCINSTALL}/lib/tgt_${G4SYSTEM}"
endif
if ( "`echo ${SHLIBVAR} | grep ${CLHEP_BASE_DIR}/lib `" == "" ) then
  set SHLIBVAR="${SHLIBVAR}:${CLHEP_BASE_DIR}/lib"
endif

setenv $SHLIBVARNAME $SHLIBVAR

#
# Remove unneeded variables. If this is not done the vars. remain in the env.
#

unsetenv SYSTEM
unset SHLIBVAR
unset SHLIBVARNAME
