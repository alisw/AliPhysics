#
# marian.ivanov@cern.ch
#
# Create  lists of missing runs
# Expected the data follow given Workspace structure
# input :  run.list  - list of runs of interest
# output:  rawMissing.list
#          esdMissing.list           
#          ocdbMissing.list
# Notice  -OCDB missing is indicated only for GRP
#         -in next version - query from the alien     

rm rawMissing.list
rm esdMissing.list
rm ocdbMissing.list 
rm rawPresent.list
rm esdPresent.list
rm ocdPresent.list 

for adir in `cat run.list`; do
   nfiles=`cat   raw$adir.txt | grep -c .root`
   if [ $nfiles -lt 2 ] ; then
    #echo RAW: 0000$adir $nfiles
    echo 0000$adir   >>rawMissing.list
    else
    echo 0000$adir   >>rawPresent.list
   fi
   nfilesReco=`cat   esd$adir.txt | grep -c .root`
   if [ $nfilesReco -lt 2 ] ; then
    #echo ESD: $adir $nfilesReco
    echo $adir   >>esdMissing.list
    else
    echo $adir   >>esdPresent.list
   fi
   nfilesOCDB=`cat grp.list | grep $adir| grep -c root`
   if [ $nfilesOCDB -lt 1 ] ; then
    #echo OCDB: $adir $nfilesOCDB
    echo $adir   >>ocdbMissing.list
    else
    echo $adir   >>ocdbPresent.list
   fi
done; 
