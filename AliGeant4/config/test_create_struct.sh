# $Id$
# ----------------------------------------------------------------
# This script generates test macro for a specified
# module and its version number with specified
# (or default) test event generator in macro/STRUCT/test/MODULE
#
# In order to prevent from unwanted rewriting of the generator
# it is created only in case it does not yet exist.
#
# Usage: 
# test_create_struct.sh modName [-d modVersionNumber] [-g genNumber] [-v visNumber]
#
# 5.11.99,  I.Hrivnacova

# check input parameters
if [ $# -lt 1 ]; then
  echo "Usage: "
  echo "test_create_struct.sh modName [-d modVersionNumber] [-g genNumber] [-v visNumber]"
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

SRC=$ALICE_ROOT/STRUCT
TO=$AG4_INSTALL/"test"/STRUCT

# test if corresponding version order number
# is defined
cd $SRC
IS_VERSION="NO"
if [ -f "Ali"$MOD"v"$VER".cxx" ]; then
  IS_VERSION="YES"
else
  if [ -f "Ali"$MOD".cxx" && "$VER" = "0" ]; then
    IS_VERSION="YES" 
  fi  
fi
if [ "$IS_VERSION" = "NO" ]; then
  cd $CURDIR
  exit;
fi       
cd $TO  

# create basic test macro
if [ "$VIS" = "0" ]; then 
  cat $AG4_INSTALL/config/test_default_det_novis.in | sed s/NNN/$VER/g | sed s/WWW/$VER/g | sed s/GGG/$GEN/g |  sed s/VVV/$VIS/g | sed s/XXX/$MOD/g > $TO/$MOD"/v"$VER"_test"$GEN""$VIS".in"
else
  cat $AG4_INSTALL/config/test_default_det_vis.in | sed s/NNN/$VER/g | sed s/WWW/$VER/g | sed s/GGG/$GEN/g |  sed s/VVV/$VIS/g | sed s/XXX/$MOD/g > $TO/$MOD"/v"$VER"_test"$GEN""$VIS".in"
# create visualisation macro (if it does not yet exist)
  if [ ! -f $TO/$MOD/vis_test$VIS".in" ]; then
    cp $AG4_INSTALL/config/test_default_vis$VIS.in $TO/$MOD/vis_test$VIS.in
  fi
fi  

# create generator macro (if it does not yet exist)
if [ ! -f $TO/$MOD/gen_test$GEN".in" ]; then
  cp $AG4_INSTALL/config/test_default_gen$GEN.in $TO/$MOD/gen_test$GEN.in
fi

cd $CURDIR
