#!/bin/bash
export ALICE_OLD=/afs/cern.ch/alice/library
#AliEn
export ALIEN=$ALICE_OLD/alien/pro
export GSHELL_ROOT=$ALIEN/api
export X509_CERT_DIR=$ALIEN/globus/share/certificates

#ROOT
source /afs/cern.ch/user/a/alisoft/.alien/packages/pchrist/ROOT/v5.15.04/.alienEnvironment /afs/cern.ch/user/a/alisoft/.alien/packages/pchrist/ROOT/v5.15.04

#AliRoot
source /afs/cern.ch/user/a/alisoft/.alien/packages/pchrist/AliRoot/v4-04-Release-a/.alienEnvironment  /afs/cern.ch/user/a/alisoft/.alien/packages/pchrist/AliRoot/v4-04-Release-a

export PATH=$PATH:$ALIEN/api/bin:$ALIEN/globus/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ALIEN/api/lib:$ALIEN/globus/lib
