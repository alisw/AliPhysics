#
# setup aliroot environment
# to be modified by users
# 
#

#
# set your aliroot
#
source /u/miranov/.balice64HEAD0108
#
#output directory
#
export GUI_OUTDIR=/lustre/alice/TPCgui

#usually the next two can stay as they are. If your not happy with where the output is written
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

#number of files per chunk in automatic tree creation
export NFILES=25

