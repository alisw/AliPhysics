# $Id$
# ----------------------------------------------------------------
# This script generates test macro for a specified
# detector and its version number with specified
# (or default) test event generator in macro/DDD/test
#
# In order to prevent from unwanted rewriting of the generator
# it is created only in case it does not yet exist.
#
# Usage: 
# create_test_in.sh detName [-d detVersionNumber] [-g genNumber] [-v visNumber]
#
# 5.11.99,  I.Hrivnacova

# check input parameters
if [ $# -lt 1 ]; then
  echo "Usage: "
  echo "create_test_in.sh detName [-d detVersionNumber] [-g genNumber] [-v visNumber]"
  exit
fi

# default parameters
DIR=$1
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

SRC=$ALICE_ROOT
TO=$AG4_INSTALL/"test"

# test if corresponding version order number
# is defined
cd $SRC/$DIR
IS_VERSION="NO"
if [ -f "Ali"$DIR"v"$VER".cxx" ]; then
  IS_VERSION="YES"
else
  if [ -f "Ali"$DIR".cxx" && "$VER" = "0" ]; then
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
  cat $AG4_INSTALL/config/test_default_det_novis.in | sed s/NNN/$VER/g | sed s/WWW/$VER/g | sed s/GGG/$GEN/g |  sed s/VVV/$VIS/g | sed s/XXX/$DIR/g > $TO/$DIR"/v"$VER"_test"$GEN""$VIS".in"
else
  cat $AG4_INSTALL/config/test_default_det_vis.in | sed s/NNN/$VER/g | sed s/WWW/$VER/g | sed s/GGG/$GEN/g |  sed s/VVV/$VIS/g | sed s/XXX/$DIR/g > $TO/$DIR"/v"$VER"_test"$GEN""$VIS".in"
# create visualisation macro (if it does not yet exist)
  if [ ! -f $TO/$DIR/vis_test$VIS".in" ]; then
    cp $AG4_INSTALL/config/test_default_vis$VIS.in $TO/$DIR/vis_test$VIS.in
  fi
fi  

# create generator macro (if it does not yet exist)
if [ ! -f $TO/$DIR/gen_test$GEN".in" ]; then
  cp $AG4_INSTALL/config/test_default_gen$GEN.in $TO/$DIR/gen_test$GEN.in
fi

cd $CURDIR
