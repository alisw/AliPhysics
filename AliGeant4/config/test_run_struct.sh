# $Id$
# ----------------------------------------------------------------
# This script runs the specified (or default) test macro 
# for a specified STRUCT module and its version number
# in aliroot/STRUCT/test/MODULE
#
# Usage: 
# test_run_struct.sh modName [-d detVersionNumber] [-g genNumber] [-v visNumber]
#
# 5.11.99,  I.Hrivnacova

# check input parameters
if [ $# -lt 1 ]; then
  echo "Usage: "
  echo "test_run_struct.sh modName [-d detVersionNumber] [-g genNumber] [-v visNumber]"
  exit
fi

# default parameters
MOD=$1
VER="0"
GEN="0"
VIS="0"

# get input parameters
for param in $*
do
  case $param in
    -d) VER=$3; shift 2;;
    -g) GEN=$3; shift 2;;
    -v) VIS=$3; shift 2;;
  esac
done

CURDIR=`pwd`
TO=$AG4_INSTALL/"test"

# go to detector test directory
cd $TO/STRUCT/$MOD

# remove old output files if exist
if [ -f "v"$VER"_test"$GEN""$VIS".out" ]; then
  rm "v"$VER"_test"$GEN""$VIS".out"
fi  
if [ -f "v"$VER"_test"$GEN""$VIS".err" ]; then
  rm "v"$VER"_test"$GEN""$VIS".err"
fi  

# run aligeant4
aligeant4 "v"$VER"_test"$GEN""$VIS".in" > "v"$VER"_test"$GEN""$VIS".out" 2> "v"$VER"_test"$GEN""$VIS".err"

# check if aligeant4 ran successfully
if [ $? -ne 0 ]; then
  echo "   !!! ERROR: Test v"$VER"_test"$GEN""$VIS" in "$MOD 
fi  
# rename g4.prim if it was created
if [ -f g4.prim ]; then
  mv g4.prim "v"$VER"_test"$GEN""$VIS".prim"
fi

cd $CURDIR
