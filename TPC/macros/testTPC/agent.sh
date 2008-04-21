#!/usr/local/bin/bash
#aguments
#1 working directory


# Example usage:
# $ALICE_ROOT/TPC/macros/testTPC/agent.sh /afs/cern.ch/user/m/miranov/public/test2008/reckr191xx

# Or submit as many agents you want on batch machines:
# bsub -q 1nd $ALICE_ROOT/TPC/macros/testTPC/agent.sh /afs/cern.ch/user/m/miranov/public/test2008/reckr191xx

#
#  AGENT to run some action:
#  3 components to be present
#

#  0. CREATE ACTION LIST ( see $ALICE_ROOT/TPC/macros/testTPC/AlienToolkit.cxx)
#

#  1. WORKING DIRECTORY SETUP (in prev. example /afs/cern.ch/user/m/miranov/public/test2008/reckr191xx)
#    1.a. action.list
#    1.b  empty out directory
#    1.c  macros directory -with  <action>.C  - e.g rec.C  FindKrClustersRaw.C


# 2. ENVIRONMENT SETUP:
#
# 2.a  Copy the script  to setup aliroot root to the  ~/.bagentsetup
# 2.b  Copy certificates to default place -  ~/.agentauth/
#     cp /tmp/*$UID* ~/.agentauth/
# 2.c  Set envirnment variables needed OCDB_PATH and INPUTTYPE


#
#     CERN SETUP - THIS WILL BE DONE IN USER ENV in future
#
export     AGENTINPUTTYPE=2 
#use       0-XRD  1-ALIEN  2-CASTOR
export     OCDB_PATH=local:///afs/cern.ch/alice/tpctest/OCDB     
#



#
#setup environment
#
#ALIEN -get token
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo ENVIRONMENT
echo ALIEN       - $ALIEN
echo GSHELL_ROOT - $GSHELL_ROOT
echo HOME        - $HOME
echo ALIROOT   =`which aliroot` 
echo    ROOT   =`which root` 
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo xxxxxxxxxxxxxxxxxxxxxxxxx
cp ~/.agentauth/*   /tmp/
source /tmp/gclient_env_$UID
export  PATH=.:${GSHELL_ROOT}/bin:${PATH}.

#aliensh -c ps
#ALIROOT
source ~/.bagentsetup
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo ENVIRONMENT
echo ALIROOT   =`which aliroot` 
echo    ROOT   =`which root` 
echo xxxxxxxxxxxxxxxxxxxxxxxxx

#
#
#
echo OPERATING SYSTEM
uname -a
#export workdir
export jobhome=`pwd`

cd $1
aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/AliTPCjobs.cxx

#alien_cp -d job.list alien:${alien_HOME}job.list@ALICE::GSI::SE 

