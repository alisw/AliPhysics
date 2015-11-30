WORKDIR=/afs/cern.ch/work/a/aliqaevs/www/data/2015/
LOG=$WORKDIR/cron.log
RUNLIST=$WORKDIR/runlist.txt
cd $WORKDIR
source env_aliroot.sh
rm $LOG
rm $RUNLIST
date > $LOG
for s in `ls /cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/GRP/GRP/LHCData/`; do
  RUN=${s:3:6}
  if (( $RUN )) 2>/dev/null; then echo $RUN >> $RUNLIST; fi
done
root -b -q runNew.C >> $LOG
aliroot -b -q integrated_lumi.C >> $LOG

