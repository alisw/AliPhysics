#
# Run train commands at GSI 
# marian.ivanov@cern.ch
#

#
# Make workspace
# This is just example (see $ALICE_ROOT/TPC/CalibMacros/ReadmeTrain.txt)
#
cp $ALICE_ROOT/TPC/macros/CalibrateTPC.C .
cp $ALICE_ROOT/TPC/macros/ConfigOCDB.C .
#modify ConfigOCDB.C
ln -sf ~/.balice64HEAD0108 balice.sh
#use your favourite aliroot
ln -sf $ALICE_ROOT/TPC/CalibMacros/alienSetupGSI.sh alienSetup.sh
#use your alien setup 
cp $ALICE_ROOT/TPC/CalibMacros/submitcalib.sh .
#copy your lists
#run.list      - main source of the run numbers  (created from the logbook)
#esdall.list   - the list of files -standard 
#esdalien.list - the list of files -alien

# Make a esd all list 
# This is just example
#
for a in `find  /lustre/alice/alien/alice/data/2009/LHC09c/0000* | grep .zip` ; do echo $a#AliESDs.root; done >esdall.list
for a in `find  /lustre/alice/alien/alice/data/2009/LHC09c/0000* | grep AliESDfriend` ; do echo `dirname $a`/AliESDs.root; done >>esdall.list



# 
# Create a list for each run 
# and make directory structure
# This is fast
mydir=`pwd`
for adir in `cat run.list`; do
    mkdir $adir;
    cat  $mydir/esdall.list | grep $adir >esd${adir}.txt
    cat  $mydir/esdalien.list | grep $adir >>esd${adir}.txt
    cp esd${adir}.txt   $adir/esd.txt
done

#
# filter lists of files 
# make 2 lists - esd.txt.Good and esd.txt.Bad
# Wait until jobs will finish
#
mydir=`pwd`
for adir in `cat run.list`; do
myvar=0;
cd $mydir
cd $adir
bsub -q alice-t3_8h  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterESD.C
echo _____________________________________
echo Run $adir
echo _____________________________________
cd $mydir
done


# You have to wait until the lists are filtered
# Only Afterwards you should process with calibration submission
# Submitting calibration
#
mydir=`pwd`
#for adir in `cat runMag5.list`; do
for adir in `cat run.list`; do
myvar=0;
cd $mydir
mkdir $mydir/$adir
cd $adir
echo SUBMITING DIRECTORY $adir
cp ../ConfigOCDB.C .
cp ../CalibrateTPC.C .
rm -rf *_*
rm -rf */V3/
up=`cat  esd.txt | grep -c .root`
while [ $myvar -le ${up} ] ; 
     do
      bsub -q alice-t3_8h  ../submitcalib.sh $myvar $(($myvar+4))  `pwd`/esd.txt.Good $adir; 
      myvar=$(( $myvar +5 )) ; 
      echo $myvar ; 
      done;
done;


#
# Submitting merging
#

mydir=`pwd`
for adir in `cat run.list`; do
myvar=0;
cd $mydir
cd $adir
rm CalibObjects*.root
ls `pwd`/*_*/V3/*/*/CalibObjectsTrain1.root > mergelistTrain1.txt
ls `pwd`/*_*/V3/*/*/CalibObjectsTrain2.root > mergelistTrain2.txt
#cat ../getmean.C | sed s_runxxx_\{$adir}\_g | sed s_{__ | sed s_}__> getmean.C
#cat ../getmeanT.C | sed "s|runxxx|$adir|" | sed s_{__ | sed s_}__> getmean.C
bsub -q alice-t3_8h aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergelistTrain1.txt\"\)
bsub -q alice-t3_8h aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergelistTrain2.txt\"\)
echo _____________________________________
echo Run $adir
echo _____________________________________
cd $mydir
done
 
#
# Merge mag field data
#
mydir=`pwd`
rm mergeMagTrain1.txt
rm mergeMagTrain2.txt
for adir in `cat runMag.list`; do
ls $mydir/$adir/CalibObjectsTrain1.root >> mergeMagTrain1.txt;
ls $mydir/$adir/CalibObjectsTrain2.root >> mergeMagTrain2.txt;
done;
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergeMagTrain1.txt\"\)
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergeMagTrain2.txt\"\)
#


#
# get debug streamer list for different components
#

find `pwd`/*/*_*/V3/ | grep calibTrigg  > trigger.txt
find `pwd`/*/*_*/V3/ | grep calibTime   > time.txt


#
# Submitting analysis
#

mydir=`pwd`
for adir in `cat esd.list`; do
myvar=0;
cd $mydir
cd $adir
ls `pwd`/*_*/V3/*/*/CalibOb* > mergelist.txt
#cat ../getmean.C | sed s_runxxx_\{$adir}\_g | sed s_{__ | sed s_}__> getmean.C
cat ../getmeanT.C | sed "s|runxxx|$adir|" | sed s_{__ | sed s_}__> getmean.C
echo _____________________________________
echo Run $adir
ls -al `pwd`/mergelist.txt
echo _____________________________________
aliroot   < getmean.C &
cd $mydir
done


ls   `pwd`/*/laserMe* > laserScan.txt 
for a in `cat laserScan.txt` ; do echo /lustre/alice/miranov/rec/LHC09d_TPC/vscan_laser1807/77475/laserMean.root; done > laserScanRefA.txt 
for a in `cat laserScan.txt` ; do echo /lustre/alice/miranov/rec/LHC09d_TPC/vscan_laser1807/77656/laserMean.root; done > laserScanRefC.txt 












