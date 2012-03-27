# Script for testing of the Pass0.jdl  on the grid for user purposes
# To be used for validation of modified user code before asking for porting

# 
# Parameters for macro:
# 1  -   alien user dir name                 - e.g  /m/miranov/
# 2  -   input directory with raw data       - e.g /alice/data/2011/LHC11a/
# 3  -   run number                          - e.g run number 000146807
# Example:
# Run=000146807
# AlienName=/m/miranov/
# RawPath=/alice/data/2011/LHC11a/
# $ALICE_ROOT/PWGPP/CalibMacros/Pass0/test/runPass0User.sh  $AlienName $RawPath  $Run 
#
# authors:   marian.ivanov#cern.ch, mikolaj.krzewicki@cern.ch 

AlienName=$1
RawPath=$2
Run=$3
InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/Pass0/

echo xxxxxxxxxxxxxxxxxxxxxxxxxx
echo SETUP
echo AlienName=$1
echo RawPath=$2
echo Run=$3 
echo InputMacros=$ALICE_ROOT/PWGPP/CalibMacros/Pass0/
echo xxxxxxxxxxxxxxxxxxxxxxxxxx

#
# 1. copy macroses and sh to the predefiend alien directory
#
OutputMacros=`echo /alice/cern.ch/user/j/jotwinow/Pass0/CalibMacros/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
alien_mkdir  $OutputMacros

for lfile in `ls $InputMacros/{*C,*sh} `; do
    bname=`basename $lfile`  
    echo  Copping alien_cp -n $lfile alien://$OutputMacros/$bname 
    alien_cp -n $lfile alien://$OutputMacros/$bname
done


#
# 2. Copy shell script and jdl
#
OutputBin=`echo  /alice/cern.ch/user/j/jotwinow/bin/ | sed s_\/j\/jotwinow\/_$AlienName\_ `
echo alien_cp -n $InputMacros/runPass0.sh  alien://$OutputBin/runPass0.sh
alien_cp -n  $InputMacros/runPass0.sh   alien://$OutputBin/runPass0.sh
cat $InputMacros/Pass0.jdl | sed "s_/j/jotwinow/_${AlienName}_g" | sed "s_/alice/data/2010/LHC10d/_${RawPath}_g" > Pass0.jdl
echo alien_cp -n Pass0.jdl alien://$OutputMacros/Pass0.jdl
alien_cp -n Pass0.jdl alien://$OutputMacros/Pass0.jdl

#
# 3. Copy validation switch off return value - job will alway finish
#
cat $InputMacros/validation.sh |  sed "s_exit \$error_exit 0_" > validation.sh
echo alien_cp  -n validation.sh  alien:///$OutputMacros/validation.sh
alien_cp  -n validation.sh  alien:///$OutputMacros/validation.sh
#
# 4. Submit job
#
echo alien_submit alien:///$OutputMacros/Pass0.jdl $Run 
#alien_submit alien:///$OutputMacros/Pass0.jdl $Run 
#echo Alien job submitted $!

