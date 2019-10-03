#!/bin/sh
# this one runs on the nodes
# works only with the correct submit script!
#
# it expects a shell var $datalistfile with the name of the file holding
# the list of files to process to be exported by the submit script
# also $listoffiles must be present.

#following makes the directories and saves the file with the files for each job
basedir=$PBS_O_WORKDIR/outputStbcAnalysis
[[ ! -d ${basedir} ]] && mkdir ${basedir}
dir=${basedir}/${PBS_JOBID}
[[ ! -d ${dir} ]] && mkdir ${dir}
echo $listoffiles
[[ -f ${dir}/${datalistfile} ]] && rm ${dir}/${datalistfile}
[[ ! -f ${dir}/${datalistfile} ]] && touch ${dir}/${datalistfile}
for x in $listoffiles
do
  echo ${x} >> ${dir}/${datalistfile}
done

#change to the proper directory ${dir}
cd ${dir}

################################################################################
## YOUR CODE HERE ##, you have $datafile in pwd to work with
################################################################################
#example:
source /project/alice/alisoft/scripts/setAlice.sh -s /data/alice3/mikolaj/alisoft/releases/trunk/
cp $PBS_O_WORKDIR/runStarFlowAnalysis.C .
root -b -q runStarFlowAnalysis.C\("\"$datalistfile\""\)

################################################################################
