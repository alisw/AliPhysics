#!/bin/sh
#ROOT
export ROOTSYS=/home/perthi/cern/root/root-current
#export ROOTSYS=/home/perthi/cern/root-current
#export ROOTSYS=/home/perthi/cern/root/root_v5.11.06
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export MANPATH=$MANPATH:$ROOTSYS/man 

#PYTHIA
export PYTHIA6=/home/perthi/cern/root/root-current/pythia6

#AliRoot
#export ALICE=/home/perthi/cern/AliRoot_head

export ALICE=/home/perthi/cern/aliroot
export ALICE_LEVEL=aliroot-current
#export ALICE=/home/perthi
#export ALICE_LEVEL=aliroot-current

export ALICE_ROOT=$ALICE/$ALICE_LEVEL
export ALICE_TARGET=`root-config --arch`
export LD_LIBRARY_PATH=$ALICE_ROOT/lib/tgt_${ALICE_TARGET}:$LD_LIBRARY_PATH
export PATH=$PATH:$ALICE_ROOT/bin/tgt_${ALICE_TARGET}


#Geant3
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH\$ALICE/geant3-current/lib/tgt_${ALICE_TARGET}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH\:$ALICE/geant3/lib/tgt_${ALICE_TARGET}
#Alice HLT components
#export ALIHLT_TOPDIR=/home/perthi/AliHLT/alihlt-current
export ALIHLT_TOPDIR=$HOME/cern/aliroot/aliroot-current/HLT
#export ALIHLT_TOPDIR=/home/perthi/cern/aliroot/aliroot-current/HLT
export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/build/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH\:$ALIHLT_LIBDIR

#HLT framework
export ALIHLT_DC_DIR=/home/perthi/HLT/hlt-current
export PATH=$PATH:$ALIHLT_DC_DIR/bin/Linux-i686
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ALIHLT_DC_DIR/lib/Linux-i686

#Xerces
export XERCESCDIR=/home/perthi/HLT/xerces
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH\:$XERCESCDIR/lib

#QT
export PLATFORM=`root-config --arch`
export QTDIR=/usr/lib/qt-3.3
export PATH=$QTDIR/bin:$PATH
export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH
export MANPATH=$QTDIR/doc/man:$MANPATH
##export LD_LIBRARY_PATH=/home/perthi/cern/aliroot/v4-04-Rev-07/HLT/build/lib:$LD_LIBRARY_PATH
##export LD_LIBRARY_PATH=/home/perthi/cern/aliroot/v4-04-Release/lib/tgt_linux:$LD_LIBRARY_PATH
#BASEDIR=/home/perthi/HLT/hlt-current
#export PATH=$PATH:$BASEDIR/bin/`uname -s`-`uname -m`
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BASEDIR/lib/`uname -s`-`uname -m`
. /home/perthi/HLT/hlt-current/bin/setenv.sh
. /home/perthi/DCS/dcs.sh
