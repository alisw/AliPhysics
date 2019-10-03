WORKDIR=/eos/user/a/aliqaevs/www/data/2017/
LOG=$WORKDIR/cron.log
RUNLIST=$WORKDIR/runlist.txt
cd $WORKDIR
source $WORKDIR/env_aliroot.sh
rm $LOG
rm $RUNLIST
date > $LOG
for s in `ls /cvmfs/alice-ocdb.cern.ch/calibration/data/2017/OCDB/GRP/GRP/LHCData/`; do
  RUN=${s:3:6}
  if (( $RUN )) 2>/dev/null; then echo $RUN >> $RUNLIST; fi
done
root -b -q runNew.C >> $LOG
