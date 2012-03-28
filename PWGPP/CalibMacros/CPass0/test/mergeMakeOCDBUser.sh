# Script for testing of the mergemakeOCDB on the grid for user purposes
# To be used in case the standard production failed or was not automatically validated
# Output objects are writn to the predefined path:

# 
# Parameters for macro:
# 1  -   alien user dir name                 - e.g  /m/miranov/
# 2  -   input directory with data to merge  - e.g /alice/data/2011/LHC11h/000169926/cpass0_HLT/
# 3  -   run number                          - e.g run number 169926
# 4  -   OCDB output path                    - alien://folder=/alice/cern.ch/user/m/miranov/Pass0/169926
# Example:
# run=170572
# AlienName=/m/miranov/
# $ALICE_ROOT/PWGPP/CalibMacros/Pass0/test/mergeMakeOCDBUser.sh $AlienName /alice/data/2011/LHC11h/000$Run/cpass0_HLT/  $Run alien://folder=/alice/cern.ch/user/m/miranov/Pass0/$Run
#
# authors:   marian.ivanov#cern.ch, mikolaj.krzewicki@cern.ch 

AlienName=$1
InputDataDir=$2
Run=$3
OCDBPath=$4
InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/Pass0/

echo xxxxxxxxxxxxxxxxxxxxxxxxxx
echo SETUP
echo AlienName=$1
echo InputDataDir=$2
echo Run=$3
echo OCDBPath=$4
echo InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/Pass0/
echo xxxxxxxxxxxxxxxxxxxxxxxxxx

#
# 1. copy macroses and sh to the predefiend alien directory
#
OutputMacros=`echo /alice/cern.ch/user/j/jotwinow/Pass0/MergeCalibration/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
alien_mkdir  $OutputMacros

for lfile in `ls $InputMacros/{*C,*sh,*jdl} | grep -v AddTask`; do
    bname=`basename $lfile`  
    echo  Copping alien_cp -n $lfile alien://$OutputMacros/$bname 
    alien_cp -n $lfile alien://$OutputMacros/$bname
done
#
# 2. Copy shell script and jdl
#
OutputBin=`echo  /alice/cern.ch/user/j/jotwinow/bin/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
echo alien_cp -n $InputMacros/mergeMakeOCDB.sh  alien://$OutputBin/mergeMakeOCDB.sh
alien_cp -n  $InputMacros/mergeMakeOCDB.sh  alien://$OutputBin/mergeMakeOCDB.sh
cat $InputMacros/mergeMakeOCDB.jdl | sed "s_/j/jotwinow/_${AlienName}_g" > mergeMakeOCDB.jdl
echo alien_cp -n mergeMakeOCDB.jdl alien://$OutputMacros/mergeMakeOCDB.jdl
alien_cp -n mergeMakeOCDB.jdl alien://$OutputMacros/mergeMakeOCDB.jdl

#
# 3. Copy validation switch off return value - job will alway finish
#
cat $InputMacros/validationMerging.sh |  sed "s_exit \$error_exit 0_" > validationMerging.sh
echo alien_cp  -n validationMerging.sh  alien:///$OutputMacros/validationMerging.sh
alien_cp  -n validationMerging.sh  alien:///$OutputMacros/validationMerging.sh
#
# 4. Submit job
#
echo alien_submit alien:///$OutputMacros/mergeMakeOCDB.jdl $InputDataDir $Run $OCDBPath
alien_submit alien:///$OutputMacros/mergeMakeOCDB.jdl $InputDataDir $Run $OCDBPath &
echo Alien job submitted $!

