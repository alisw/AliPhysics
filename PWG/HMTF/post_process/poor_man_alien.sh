#!/bin/bash

# export LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/":$LD_LIBRARY_PATH
# export PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin/":$PATH 
export LD_LIBRARY_PATH="/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliEn-Runtime/v2-19-le-18/lib/":$LD_LIBRARY_PATH
export PATH="/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliEn-Runtime/v2-19-le-18/bin/":$PATH
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6-orig/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.30/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
# alien_cp checks for gbbox by using GSHELL_ROOT...
export GSHELL_ROOT="/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliEn-Runtime/v2-19-le-18"
export X509_CERT_DIR=$GSHELL_ROOT"/globus/share/certificates"

# source python virtual env
mkdir -p ~/.virtualenvs
virtualenv -q ~/.virtualenvs/hmtfc
source .virtualenvs/hmtfc/bin/activate
pip install --upgrade -e post_process/
