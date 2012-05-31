#
# Sequence of steps to test Pass0 and PassX reconstruction/calibration which run on GRID 
# by default               
# 
# Semi automatic test performed on the batch farm
# Important features:
# 1. Parsing of the log files                           
# 2. Parsing stack traces
   
# author:  marian.ivanov@cern.ch


#
# 0. copy a template setup
#   To be modified if necessary
#  "standard scripts
cp $ALICE_ROOT/ANALYSIS/macros/runCalibTrain.C  .
#  scipts to run  on batch farm
cp $ALICE_ROOT/TPC/scripts/calibPassX/recPass0GSI.C   rec.C 
cp $ALICE_ROOT/TPC/scripts/calibPassX/runPassX.sh     .
cp $ALICE_ROOT/TPC/scripts/calibPassX/submitPass0.sh  .
cp $ALICE_ROOT/TPC/scripts/calibPassX/runPassXJob.sh  .

cp ../lists/run.list .
cp ../lists/raw.list .

#
# 1. Make workspace - directory structure with run and raw lists 
#
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list 
#
# 1.a clean the workspace all
#
find `pwd` | grep AliESD  | xargs rm
find `pwd` | grep out     | xargs rm
find `pwd` | grep log     | xargs rm
#
# 1.b clean workscpace  - rm root files if not tags found
#
find `pwd` | grep AliESDs.root > esd.txt
for efile in `cat esd.txt`; do
    dname=`dirname $efile` 
    status=`ls $dname/*tag.root`; 
    echo  CHECK  $efile $status
    if [ -z $status ]; then 
      echo NON OK  rm -r $dname/*.root
      rm -r $dname/*.root
    else
      echo IS OK $status
    fi; 
done

#
# 2. Run reconstruction/calibration
#
bgroup=/recPass0/`pwd | xargs basename`
bgadd $bgroup
nEvents=500000
submitPass0.sh run.list "alice-t3_8h -g  $bgroup" $nEvents | tee  submit.log


#
# 3. Run merging - run level
#
for arun in `cat run.list`; do
    cd $arun
    find `pwd`/ | grep AliESDfriends_v1.root > calib.list
    echo bsub -q alice-t3_8h  -o outMerge.log  aliroot -b -q  $ALICE_ROOT/ANALYSIS/macros/mergeCalibObjects.C
    bsub -q alice-t3_8h -o outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/macros/mergeCalibObjects.C
    cd ../
done;

#
# 4. Merge all
#
ls `pwd`/*/CalibObjects.root  > calib.list
aliroot -b -q  $ALICE_ROOT/ANALYSIS/macros/mergeCalibObjects.C

#
# 5. filter reconstruction  logs, make statistic of problems
#
mkdirhier logs/reco
cd  logs/reco
find `pwd`/../../ | grep rec.log > errRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh
cd ../..

#
# 6. filter calibration  logs, make statistic of problems
#
mkdirhier logs/calib
cd  logs/calib
find `pwd`/../../ | grep calib.log > errRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh
cd ../..


#
# 7. filter debug streamer
#
mkdir debug
cd debug
find  `pwd`/../*/ | grep V6 | grep .root  > debugall.txt
cat  debugall.txt | grep calibTimeDebug > timeitstpc.txt




