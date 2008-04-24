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
#export     AGENTINPUTTYPE=2 
#use       0-XRD  1-ALIEN  2-CASTOR
#export     OCDB_PATH=local:///afs/cern.ch/alice/tpctest/OCDB     
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

#try to cd to workdir
res=`cd $1 || echo no`
if [ "x$res" == "xno" ]; then
  echo "could not cd to workdir '$1'"
  exit 1
fi

# check if we are able to write to this dir, the alien SE, can readback from the alien SE and delete stuff
touch agend_alien_tst.txt || exit 1
date >> agend_alien_tst.txt
#copy to alien
alien_cp -d agend_alien_tst.txt alien:${alien_HOME}agend_alien_tst.txt@$AGENTSE

#see if it is there
res=`alien_ls ${alien_HOME}agend_alien_tst.txt | grep agend_alien_tst`
if [ "x$res" == "x" ]; then
  echo "could not write on storage element '$AGENTSE'"
  exit 1
fi

#copy back
alien_cp -n alien:${alien_HOME}agend_alien_tst.txt agend_alien_tst_back.txt 
res=`diff agend_alien_tst.txt agend_alien_tst_back.txt`
if [ "x$res" != "x" ]; then
  echo "problems reading from storage element '$AGENTSE'"
  exit 1
fi

#try to delete from alien
res=`alien_rm -d ${alien_HOME}agend_alien_tst.txt`
if [ "x$res" != "x" ]; then
  echo "cannot delete from storage element '$AGENTSE'"
  exit 1
fi

#try to delete local files
res=`rm agend_alien_tst*`
if [ "x$res" != "x" ]; then
  echo "cannot delete local files"
  exit 1
fi



aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/AliTPCjobs.cxx

#alien_cp -d job.list alien:${alien_HOME}job.list@$AGENTSE

