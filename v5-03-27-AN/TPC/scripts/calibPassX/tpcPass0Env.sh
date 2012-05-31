#
# parameters:
# 1 - basedir
# Example: 
# source  /usr/local/grid/AliRoot/HEAD0108/TPC/scripts/tpcPass0Env.sh `pwd`
export balice=/u/miranov/.balice
export aliensetup=$HOME/alienSetup.sh
export PASS0_DIR=/usr/local/grid/AliRoot/HEAD0108
source $balice
#source $aliensetup >aliensetup.log
#
# Test setup
#
export workdir=$1
if [ ! -n length ]; then 
  echo \############################  
  echo Directory was not specified. Exiting
  echo \############################   
  return;
fi;
if [ ! -r $workdir/lists/esd.list  ] ; then
 echo \############################   
 echo File esd list does not exist. Exiting
 echo \############################   
 return;
fi; 
if [ ! -r $workdir/lists/run.list  ] ; then
 echo \############################   
 echo File run list does not exist. Exiting
 echo \############################   
 return;
fi; 

#
# Make directories
#
cd $workdir
chgrp -R alice $workdir
chmod -R g+rwx $workdir
chmod -R o+rx $workdir
mkdirhier  $workdir/calibNoDrift
mkdirhier  $workdir/calibNoRefit
mkdirhier  $workdir/calibQA
#
#modify ConfigOCDB.C
#
# copy predefined Config files 
#
cp   $PASS0_DIR/TPC/macros/CalibrateTPC.C      calibNoDrift/CalibrateTPC.C
cat  $PASS0_DIR/TPC/macros/CalibrateTPC.C |    grep -v AddCalibCalib\(task\) > calibNoRefit/CalibrateTPC.C
cp   $PASS0_DIR/TPC/macros/CalibrateTPC.C      calibQA/CalibrateTPC.C
cp   $PASS0_DIR/TPC/macros/ConfigOCDBNoDrift.C calibNoDrift/ConfigOCDB.C
cp   $PASS0_DIR/TPC/macros/ConfigOCDBNoRefit.C calibNoRefit/ConfigOCDB.C
cp   $PASS0_DIR/TPC/macros/ConfigOCDBQA.C      calibQA/ConfigOCDB.C
cp   lists/*.list calibNoDrift/
cp   lists/*.list calibNoRefit/
cp   lists/*.list calibQA/
ln -sf $balice          calibNoDrift/balice.sh
ln -sf $balice          calibNoRefit/balice.sh
ln -sf $balice          calibQA/balice.sh
ln -sf $aliensetup      calibNoDrift/alienSetup.sh
ln -sf $aliensetup      calibNoRefit/alienSetup.sh
ln -sf $aliensetup      calibQA/alienSetup.sh
#  make workspaces
#
cd $workdir/calibNoDrift
$PASS0_DIR/TPC/scripts/makeWorkspace.sh run.list 
$PASS0_DIR/TPC/scripts/submitCalib.sh run.list alice-t3 20
cd $workdir/calibNoRefit
$PASS0_DIR/TPC/scripts/makeWorkspace.sh run.list 
$PASS0_DIR/TPC/scripts/submitCalib.sh run.list alice-t3 20
cd $workdir/calibQA
$PASS0_DIR/TPC/scripts/makeWorkspace.sh run.list 
$PASS0_DIR/TPC/scripts/submitCalib.sh run.list alice-t3 20
cd $workdir/
#
#
#
