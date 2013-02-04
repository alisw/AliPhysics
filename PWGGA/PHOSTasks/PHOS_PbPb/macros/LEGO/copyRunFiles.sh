# This macro can be used to copy the run merge files of a lego train to local machine
# You may need to change some of the following variables
# author: Henrik Qvigstad <henrik.qvigstad@cern.ch>

# From where to read:
RUNLISTFULL=runlistLHC11h.txt
AFILENAME=AnalysisResults.root
GRID_PATH=/alice/data/2011/LHC11h_2
GRID_FILE_PATERN=ESDs/pass2/AOD115/PWGGA/GA_PbPb/17_20121218-1528/$AFILENAME
# Where to store:
LOCAL_DIR=$(pwd)/GA_PbPb_17
FIND_RESULTS_FILE=$LOCAL_DIR/find_results.txt
RUNFILE=$LOCAL_DIR/runFile.txt



if [ -f "$FIND_RESULTS_FILE" ]; then
    echo using $FIND_RESULTS_FILE
else
    mkdir -p $LOCAL_DIR
    echo Doing alien_find $GRID_PATH $GRID_FILE_PATERN, 
    echo may take some time ...
    alien_find $GRID_PATH $GRID_FILE_PATERN >> $FIND_RESULTS_FILE
    echo done
    echo storing results to $FIND_RESULTS_FILE
fi

# Remove the file which lists runs and files
rm -rf $RUNFILE

# Download run output files, and fill $RUNFILE
for run in $(cat $RUNLISTFULL); do
    echo looking for run $run in find results
    for line in $(cat $FIND_RESULTS_FILE); do
	file=$(echo $line | grep $run)
	if [ -n "$file" ]; then
	    break
	fi
    done
    if [ -n "$file" ]; then
	TOPATH=$LOCAL_DIR/$run
	TOFILE=$TOPATH/$AFILENAME
	mkdir -p $TOPATH
	if [[ -f "$TOFILE" && -s "$TOFILE" ]]; then
	    echo file exists, abort copying.
	else
	    echo alien_cp alien:$file file:$TOFILE
	    alien_cp alien:$file file:$TOFILE
	fi
	echo adding to $RUNFILE 
	echo "    $run and $TOFILE"
	echo "$run $TOFILE" >> $RUNFILE
    else
	echo $run has no file
    fi
done
