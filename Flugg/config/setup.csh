# $Id$
# Flugg tag $Name$
# ----------------------------------------------------------------
# This script sets the default environment variables for
# ALICE Flugg;
# the Geant4 environment is supposed to be set by setup.(c)sh
# script provided in AliRoot.  
#   
# by I. Hrivnacova, 10.5.2001
# csh version by I.Gonzalez  6.3.2002 

# check ALICE Geant4 environment
if ( "$AG4_INSTALL" == ""  ) then
  echo "ALICE Geant4 environment is not defined."
  echo "Check your AliRoot environment and run:"
  echo "source $ALICE_ROOT/AliGeant4/config/setup.csh"
  exit 1;
endif

# FLUKA environment
setenv FLUKA /afs/cern.ch/alice/offline/fluka/01c
setenv FORLIB $FLUKA/libflukahp.a
setenv FLUPRO $FLUKA

# Flugg environment
setenv FLUGGINSTALL ${ALICE_ROOT}/Flugg
#setenv G4GEOMETRY_DEBUG=1

# path to Flugg shared libararies
if ( "$G4LIB_BUILD_SHARED" != "" ) then
  set SYSTEM=`uname`
  if ( "$SYSTEM" == "Linux" ) then
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${FLUGGINSTALL}/lib/${G4SYSTEM}"
  endif
  if ( "$SYSTEM" == "HP-UX" ) then
    setenv SHLIB_PATH "${SHLIB_PATH}:${FLUGGINSTALL}/lib/${G4SYSTEM}"
  endif
  if ( "$SYSTEM" == "OSF1" ) then
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${FLUGGINSTALL}/lib/${G4SYSTEM}"
  endif
  if ( "$SYSTEM" == "SunOS" ) then
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${FLUGGINSTALL}/lib/${G4SYSTEM}"
  endif
endif
# path to Flugg executable
setenv PATH "${PATH}:${FLUGGINSTALL}/bin/Linux-g++"

#Example of run command:
#$FLUKA/flutil/rfluka -M1 -N0 -p wa_50m 
#-e $FLUGGINSTALL/bin/Linux-g++/mainAlAuAl alaual

unset SYSTEM
