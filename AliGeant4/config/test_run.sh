# $Id$
# ----------------------------------------------------------------
# This script runs specified (or default) test macros 
# for all versions of one detector if it specified 
# or for all detectors otherwise
#
# Usage: test_run.sh detName [-g genNumber] [-v visNumber]
#                      for all detectors: detName = ALL
#
# 5.11.99  I. Hrivnacova

# check input parameters
if [ $# -lt 1 ]; then
  echo "Usage: "
  echo "test_run.sh detName [-g genNumber] [-v visNumber]"
  echo "               for all detectors: detName = ALL"
  exit
fi

# default parameters
DIR=$1
GEN="0"
VIS="0"

# get input parameters
for param in $*
do
  case $param in
    -g) GEN=$3; shift 2;;
    -v) VIS=$3; shift 2;;
  esac
done

CURDIR=`pwd`
SRC=$ALICE_ROOT
TO=$AG4_INSTALL/"test"
MAX=10

cd $TO
if [ "$DIR" = "ALL" ]; then
# loop over all detectors if det is not specified
  for DIR in `ls`; do
    if [ -d $DIR ] ; then
      cd $TO/$DIR
      VER=0
      until [ "$VER" = "$MAX" ] ; do
        if [ -f "v"$VER"_test"$GEN""$VIS".in" ] ; then
          echo "test_run_in.sh $DIR v$VER test$GEN$VIS"
          test_run_in.sh $DIR -d $VER -g $GEN -v $VIS
        fi
        let VER=$VER+1
      done 
      cd $TO
    fi       
  done
else

# loop over all structures
  if [ "$DIR" = "STRUCT" ]; then
    cd $TO/$DIR
    for MODULE in `ls`; do
      VER=0
      until [ "$VER" = "$MAX" ] ; do
        if [ -f $MODULE/"v"$VER"_test"$GEN""$VIS".in" ]; then
          echo "test_run_struct.sh $MODULE v$VER test$GEN$VIS"
          test_run_struct.sh $MODULE -d $VER -g $GEN -v $VIS
	fi  
        let VER=$VER+1
      done
    done
  else
  
# run for specified detector only
    if [ -d $DIR ] ; then
      cd $TO/$DIR
      VER=0
      until [ "$VER" = "$MAX" ] ; do
        if [ -f "v"$VER"_test"$GEN""$VIS".in" ] ; then
          echo "test_run_in.sh $DIR v$VER test$GEN$VIS"
          test_run_in.sh $DIR -d $VER -g $GEN -v $VIS
        fi
        let VER=$VER+1
      done 
      cd $TO
    fi       
  fi  
fi

cd $CURDIR
