#!/bin/sh
# $Id$
# ----------------------------------------------------------------
# This script generates test macros for all versions of
# one detector if it specified or for all detectors otherwise
# with specified (or default) test event generator 
# in macro/DDD/test
#
# Usage: test_create.sh detName [-g genNumber] [-v visNumber]
#                      for all detectors: detName = ALL
#
# 5.11.99  I. Hrivnacova

# check input parameters
if [ $# -lt 1 ]; then
  echo "Usage: "
  echo "test_create.sh detName [-g genNumber] [-v visNumber]"
  echo "                 for all detectors:  detName = ALL"
  echo "                     all structures: detName = STRUCT"
  echo "                     all detectors of all versions:  detName = ALL_VERSIONS"
  echo "                     all structures of all versions: detName = STRUCT_ALL_VERSIONS"
  exit
fi

# default parameters
DIR=$1
GEN="1"
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

# create destination directory if it does not exist
if [ ! -d $TO ] ; then    
  mkdir $TO
  cd $TO
  mkdir CRT EMCAL FMD ITS MUON PHOS PMD RICH START STRUCT TOF TPC TRD ZDC
  cd STRUCT
  mkdir ABSO DIPO FRAME HALL MAG PIPE SHIL
fi  

cd $SRC
if [ "$DIR" = "ALL_VERSIONS" ]; then
# loop over all detectors if det is not specified
  for DIR in `ls`; do
    if [ -d $DIR ] ; then
      cd $SRC/$DIR
      VER=0
      until [ "$VER" = "$MAX" ] ; do
        if [ -f "Ali"$DIR"v"$VER".cxx" ]; then
          #echo "test_create_in.sh $DIR v$VER test$GEN$VIS"
          test_create_in.sh $DIR -d $VER -g $GEN -v $VIS
        fi
        let VER=$VER+1
      done 
      cd $SRC    
    fi  
  done
else

# loop over all structures
  if [ "$DIR" = "STRUCT_ALL_VERSIONS" ]; then
    cd $TO/$DIR
    for MODULE in `ls`; do
      VER=0
      until [ "$VER" = "$MAX" ] ; do
        if [ -f $SRC/$DIR/"Ali"$MODULE"v"$VER".cxx" ]; then
          #echo "test_create_struct.sh $MODULE v$VER test$GEN$VIS"
          test_create_struct.sh $MODULE -d $VER -g $GEN -v $VIS
	fi  
        let VER=$VER+1
      done
    done
  else
  
# all detectors with default versions only
    if [ "$DIR" = "ALL" ]; then
      cd $SRC
      for DIR in `ls`; do
        if [ -d $DIR ]; then
          cd $SRC/$DIR
          VER="d"
          #echo "test_create_in.sh $DIR v$VER test$GEN$VIS"
          test_create_in.sh $DIR -d $VER -g $GEN -v $VIS
        fi
        cd $SRC    
      done
    else

# loop over all structures with default versions only
      if [ "$DIR" = "STRUCT" ]; then
        DIR="STRUCT"
        cd $TO/$DIR
        for MODULE in `ls`; do
          VER="d"
          #echo "test_create_struct.sh $MODULE v$VER test$GEN$VIS"
          test_create_struct.sh $MODULE -d $VER -g $GEN -v $VIS
        done
      else
  
# specified detector only
        if [ -d $SRC/$DIR ]; then    
          cd $SRC/$DIR
          VER=0
          until [ "$VER" = "$MAX" ] ; do
            if [ -f "Ali"$DIR"v"$VER".cxx" ]; then
              #echo "test_create_in.sh $DIR v$VER test$GEN$VIS"
              test_create_in.sh $DIR -d $VER -g $GEN -v $VIS
            fi
            let VER=$VER+1
          done 
          cd $TO     
        fi
      fi
    fi
  fi    
fi
cd $CURDIR
