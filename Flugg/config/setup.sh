# $Id$
# Flugg tag $Name$
# ----------------------------------------------------------------
# This script sets the default environment variables for
# ALICE Flugg;
# the Geant4 environment is supposed to be set by setup.(c)sh
# script provided in AliRoot.  
#   
# by I. Hrivnacova, 10.5.2001

# check ALICE Geant4 environment
if [ "$AG4_INSTALL" = ""  ]; then
  echo "ALICE Geant4 environment is not defined."
  echo "Check your AliRoot environment and run:"
  echo ". $ALICE_ROOT/AliGeant4/config/setup.sh"
  exit 1;
fi

# FLUKA environment
export FLUKA=/afs/cern.ch/alice/offline/fluka/01c
export FORLIB=$FLUKA/libflukahp.a
export FLUPRO=$FLUKA

# Flugg environment
export FLUGGINSTALL=$ALICE_ROOT/Flugg
#export G4GEOMETRY_DEBUG=1

# path to Flugg shared libararies
if [ "$G4LIB_BUILD_SHARED" != "" ]; then
  SYSTEM=`uname`
  if [ "$SYSTEM" = "Linux" ]; then
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FLUGGINSTALL/lib/$G4SYSTEM
  fi
  if [ "$SYSTEM" = "HP-UX" ]; then
    SHLIB_PATH=$SHLIB_PATH:$FLUGGINSTALL/lib/$G4SYSTEM
  fi
  if [ "$SYSTEM" = "OSF1" ]; then
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FLUGGINSTALL/lib/$G4SYSTEM
  fi
  if [ "$SYSTEM" = "SunOS" ]; then
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FLUGGINSTALL/lib/$G4SYSTEM
  fi
fi
# path to Flugg executable
export PATH=$PATH:$FLUGGINSTALL/bin/Linux=g++


#Example of run command:
#$FLUKA/flutil/rfluka -M1 -N0 -p wa_50m 
#-e $FLUGGINSTALL/bin/Linux-g++/mainAlAuAl alaual

unset SYSTEM
