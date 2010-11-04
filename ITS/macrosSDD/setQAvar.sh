#!/bin/bash

#environment script to launch alice macros on lxplus from a script of Renu Bala (Thanks Renu for the help)

export ALICE_OLD=/afs/cern.ch/alice/library

#AliEn
export ALIEN=$ALICE_OLD/alien/pro
export GSHELL_ROOT=$ALIEN/api
export X509_CERT_DIR=$ALIEN/globus/share/certificates

#ROOT
source /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/ROOT/v5-26-00b-6/.alienEnvironment /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/ROOT/v5-26-00b-6
source /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/GEANT3/v1-11-11/.alienEnvironment  /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/GEANT3/v1-11-11/

#AliRoot
source /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/AliRoot/v4-19-18-AN/.alienEnvironment  /afs/cern.ch/user/a/alisoft/.alien/packages/VO_ALICE/AliRoot/v4-19-18-AN
export PATH=$PATH:$ALIEN/api/bin:$ALIEN/globus/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ALIEN/api/lib:$ALIEN/globus/lib

alias launchalien=$GSHELL_ROOT/bin/alien-token-init
