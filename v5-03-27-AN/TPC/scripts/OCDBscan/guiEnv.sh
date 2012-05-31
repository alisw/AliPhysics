#
# setup aliroot environment
# to be modified by users
# This is example setup which is used at GSI
# In order to use it on your laptop 
# AliRoot, Destination directories and the AUTOFILES has to be modified
# (Jens Wiechula, Marian Ivanov)
#
#
# parameters
# 1  - debug flag

#
# set your aliroot and  alien environment
#
#source $HOME/.balice64HEAD0108
#source $HOME/.balice
echo YOU HAVE TO MODIFY ALIROOT SETUP
echo IT  IS ENVIRNMENT SPECIFIC 
#
#output directory - to be set according to your setup
#
export GUI_OUTDIR=/home/kowalski/alice/PPL
echo YOU HAVE TO MODIFY DESTINATION DIRECTORY
echo IT  IS ENVIRONMENT SPECIFIC 

#usually the next two can stay as they are. If you are not happy with where the output is written
#feel free to change them
export GUI_OUTDIR_TIME=$GUI_OUTDIR/time
export GUI_OUTDIR_RUNS=$GUI_OUTDIR/runs

#command for batch processing
# see also TMPLISTDIR!!!
export BATCHCOMMAND="bsub -q alice-t3_8h"

#directory for the temporary list files
#for batch processing this dir needs to be accessable from the batch nodes!!!
export TMPLISTDIR=$GUI_OUTDIR/guiTreeLists

#path to reference tree. Consider to use one!
#see README to understand how to create them
export REF_DATA_FILE=$GUI_OUTDIR/ref/RefCalPads.root
export REF_DATA_TREE=$GUI_OUTDIR/ref/RefTree.root

#whether alien path is used in OCDB
#export WITHALIEN=1
export WITHALIEN=0

#where to look for files in automatic tree creation
#it assumes a path to an OCDB directory and uses the
#run numbers of the file names:
#Run([0-9]{5})_.*
export AUTOFILES=/lustre/alice/alien/alice/data/2009/OCDB/TPC/Calib/HighVoltage
echo YOU HAVE TO MODIFY AUTOFILES DIRECTORY
echo IT  IS ENVIRONMENT SPECIFIC 

#number of files per chunk in automatic tree creation
export NFILES=25
 


echo Test guiEnv setup

errorCode=0;
goodPass=0;
testROOT=`which root`
if [ -z testROOT ];  then
   echo Check root setup"      ":FALSE
   let errorCode=errorCode+1
 else
   echo Check root setup"      ":OK"    "-  $testROOT
fi;
#
# test GUIdir
#
if [ -d $GUI_OUTDIR ];  then
   echo Check GUI_OUTDIR"      ":OK"    "-  $GUI_OUTDIR
  else
   echo Check GUI_OUTDIR"      ":FALSE"" -  $GUI_OUTDIR
   let errorCode=errorCode+2
fi;
#
#
#
if [ -r $GUI_OUTDIR_TIME ];  then
    echo Check GUI_OUTDIR_TIME" ":OK"    "-  $GUI_OUTDIR_TIME
  else
    echo Check GUI_OUTDIR_TIME" ":FALSE" "-  Does not exist or not readable           
    let errorCode=errorCode+4
fi;

if [ -r $GUI_OUTDIR_RUNS ];  then
   echo Check GUI_OUTDIR_RUNS" ":OK"    "-  $GUI_OUTDIR_RUNS
  else
   echo Check GUI_OUTDIR_RUNS" ":FALSE" "-  Does not exist or not readable           
   let errorCode=errorCode+8
fi;

if [ -r $GUI_OUTDIR/guiTreeLists ];  then
   echo Check GUI_OUTDIR/guiTreeLists" ":OK"    "-  $GUI_OUTDIR/guiTreeLists
  else
   echo Check GUI_OUTDIR/guiTreeLists" ":FALSE" "-  Does not exist or not readable           
   let errorCode=errorCode+16
fi;

if [ -z `which aliensh` -a WITHALIEN!=0 ]; then
   echo Alien not properly initialized
   let errorCode=errorCode+32
fi;


echo $errorCode
