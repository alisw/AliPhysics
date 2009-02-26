#!/bin/sh

# 1 argument      - the path to the environment setup
# 2 argument      - the job ID
# 3 argument      - path to Config.C  file
# 4 argument      - path to the database files
# 5 argument      - number of events in the file
# 6 argument      - output path
# 7 argument      - reco type
# 8 argument      - working directory


# EXAMPLE
# 1.
# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 2 $ALICE_ROOT/TPC/testMC/ConfigHPT.C  0  2 hpt 0 `pwd`

# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 2 $ALICE_ROOT/TPC/testMC/ConfigPP.C  0  10 pp 0 `pwd`

# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 2 $ALICE_ROOT/TPC/testMC/ConfigCosmic.C  0  10 cosmic 0 `pwd`


# myvar=0
# while [ $myvar -ne 100 ] ; do bsub  do something ;  myvar=$(( $myvar + 1 )) ; echo $myvar ; done

# 1.b 
# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 0  $ALICE_ROOT/TPC/testMC/ConfigHPT1.C  0  10 hpt1 0 `pwd`
# 2.
# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 0  $ALICE_ROOT/TPC/testMC/Config_AliGenCosmicsParam.C   0  2 cosmic 2
# 3. $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 0  $ALICE_ROOT/TPC/testMC/ConfigHM.C  0  1 hm 0
# $ALICE_ROOT/TPC/testMC/submitMC.sh /u/miranov/.balice64HEAD0108 0  $ALICE_ROOT/TPC/testMC/ConfigLM.C  0  1 lm 0



#
# 1  /u/miranov/.balice64v4-06-Release   # setup aliroot -root
# 2  0                                   # local directory path
# 3 \$ALICE_ROOT/macros/ConfigHPT.C      # path to the Config file
# 4 0                                    # path to particular TPC calib files
# 5 2                                    # number of events per file
# 6 hpt                                  # path where wi will write the output


cd $8
mkdir $6
cd $6
cp $3 .
cp $ALICE_ROOT/TPC/testMC/sim.C .
cp $ALICE_ROOT/TPC/testMC/recMC.C .



echo PWD `pwd` 

echo HOSTNAME $HOSTNAME
# 1 SETUP given ROOT and ALIROOT
echo   $1
source $1
echo  $ROOTSYS
which root.exe
which aliroot
  
mkdir $2
cd $2
cp ~/rootlogon.C .
echo Job ID  $2
echo
echo PWD `pwd`
 
#
#
#####################################################################
echo SUBMITING MACRO
echo "$ALICE_ROOT/TPC/testMC/sim.C(\"$3\",\"$4\",$5)"
echo
command aliroot  -q -b "$ALICE_ROOT/TPC/testMC/sim.C(\"$3\",\"$4\",$5)"
echo
echo 
#
#
#
rm *.Hits*root
rm *.SDi*root
#####################################################################
echo 
echo SUBMITING  RECONSTRUCTION MACRO
echo "$ALICE_ROOT/TPC/testMC/recMC.C"
rm AliESD*
command aliroot  -q -b  "$ALICE_ROOT/TPC/testMC/recMC.C(\"$4\")"


