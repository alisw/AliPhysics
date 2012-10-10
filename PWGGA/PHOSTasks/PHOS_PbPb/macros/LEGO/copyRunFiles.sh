# This macro can be used to copy the run merge files of a lego train to local machine
# You may need to change some of the following variables
# author: Henrik Qvigstad <henrik.qvigstad@cern.ch>

# From where to read:
RUNLISTFULL=runlistLHC11h.txt
AFILENAME=AnalysisResults.root
GRID_PATH=/alice/data/2011/LHC11h_2
GRID_FILE_PATERN=ESDs/pass2/AOD095/PWGGA/GA_PbPb/5_20121003-2356/$AFILENAME
# Where to store:
LOCAL_DIR=$(pwd)/GA_PbPb_output_5
RUNLIST=$LOCAL_DIR/runlist.txt
FILELIST=$LOCAL_DIR/filelist.txt  




echo Doing alien_find $GRID_PATH $GRID_FILE_PATERN, may take some time ...
FIND_RESULTS=$(alien_find $GRID_PATH $GRID_FILE_PATERN)
echo done, copying files.

rm -f $RUNLIST
rm -f $FILELIST

for run in $(cat $RUNLISTFULL); do
    echo looking for run $run in find results
    for line in $FIND_RESULTS; do
	file=$(echo $line | grep $run)
	if [ -n "$file" ]; then
	    break
	fi
    done
    if [ -n "$file" ]; then
	TOPATH=$LOCAL_DIR/$run
	TOFILE=$TOPATH/$AFILENAME
	mkdir -p $TOPATH
	if [ -f $TOFILE ]; then
	    echo file exists, abort copying.
	else
	    echo alien_cp alien:$file file:$TOFILE
	    alien_cp alien:$file file:$TOFILE
	fi
	echo adding to $RUNLIST $run and to $FILELIST $TOFILE
	echo $run >> $RUNLIST
	echo $TOFILE >> $FILELIST
    else
	echo $run has no file
    fi
done
