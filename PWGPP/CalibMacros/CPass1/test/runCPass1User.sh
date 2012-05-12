# Script for testing of the CPass1.jdl  on the grid for user purposes
# To be used for validation of modified user code before asking for porting

# 
# Parameters for macro:
# 1  -   alien user dir name                 - e.g  /m/miranov/
# 2  -   input directory with raw data       - e.g /alice/data/2011/LHC11a/
# 3  -   run number                          - e.g run number 146807
# Example:
# RunNumber=146807
# AlienName=/m/miranov/
# RawPath=/alice/data/2011/LHC11a/
# $ALICE_ROOT/PWGPP/CalibMacros/CPass1/test/runCPass1User.sh  $AlienName $RawPath  $RunNumber 
#
# authors:   marian.ivanov#cern.ch, mikolaj.krzewicki@cern.ch 

AlienName=$1
RawPath=$2
RunNumber=$3
InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/CPass1/

echo xxxxxxxxxxxxxxxxxxxxxxxxxx
echo SETUP
echo AlienName=$1
echo RawPath=$2
echo RunNumber=$3 
echo InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/CPass1/
echo xxxxxxxxxxxxxxxxxxxxxxxxxx

#
# 1. copy macroses and sh to the predefiend alien directory
#
OutputMacros=`echo /alice/cern.ch/user/j/jotwinow/CPass1/CalibMacros/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
alien_mkdir -p $OutputMacros

for lfile in `ls $InputMacros/{*C,*sh} `; do
    bname=`basename $lfile`  
    echo  Copping alien_cp -n $lfile alien://$OutputMacros/$bname 
    alien_cp -n $lfile alien://$OutputMacros/$bname
done


#
# 2. Copy shell script and jdl
#
OutputBin=`echo  /alice/cern.ch/user/j/jotwinow/bin/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
echo alien_cp -n $InputMacros/runCPass1.sh  alien://$OutputBin/runCPass1.sh
alien_cp -n  $InputMacros/runCPass1.sh   alien://$OutputBin/runCPass1.sh
cat $InputMacros/CPass1.jdl | sed "s_/j/jotwinow/_${AlienName}_g" | sed "s_/alice/data/2010/LHC10d/_${RawPath}_g" > CPass1.jdl
echo alien_cp -n CPass1.jdl alien://$OutputMacros/CPass1.jdl
alien_cp -n CPass1.jdl alien://$OutputMacros/CPass1.jdl

#
# 3. Copy validation switch off return value - job will alway finish
#
cat $InputMacros/validation.sh |  sed "s_exit \$error_exit 0_" > validation.sh
echo alien_cp  -n validation.sh  alien:///$OutputMacros/validation.sh
alien_cp  -n validation.sh  alien:///$OutputMacros/validation.sh
#
# 4. Submit job
#
echo nohup alien_submit alien:///$OutputMacros/CPass1.jdl "000"$RunNumber  >submitJob$RunNumber.txt
nohup alien_submit alien:///$OutputMacros/CPass1.jdl "000"$RunNumber  >submitJob$RunNumber.txt
#echo Alien job submitted $!

