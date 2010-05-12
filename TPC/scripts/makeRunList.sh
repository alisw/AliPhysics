#
# Make a run list using
# GRP and TPC HV and Altro requested
# Output run list stored  in run list

prefix=/alice/data/2010/OCDB/

alien_find $prefix/GRP/GRP/Data Run > grpAlien.txt
alien_find $prefix/TPC/Calib/HighVoltage Run > hvAlien.txt
alien_find $prefix/TPC/Calib/AltroConfig Run > altroAlien.txt

cat grpAlien.txt | sed s/_/\ /g | gawk '{ print $2}' | sort > grp.txt
cat hvAlien.txt | sed s/_/\ /g | gawk '{ print $2}'  | sort > hv.txt
cat altroAlien.txt | sed s/_/\ /g | gawk '{ print $2}'  | sort > altro.txt

for run in `cat hv.txt | sort`; do
   grun=`cat grp.txt | grep -c $run`
   arun=`cat altro.txt | grep -c $run`   
   if [ $grun -gt 0 ] && [ $arun -gt 0 ]; then
      echo $run
   fi 
done > run.list


