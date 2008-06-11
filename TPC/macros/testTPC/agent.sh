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
#    1.a  action.list
#    1.b  empty out directory
#    1.c  macros directory -with  <action>.C  - e.g rec.C  FindKrClustersRaw.C


# 2. ENVIRONMENT SETUP:
#
# 2.a  Copy the script  to setup aliroot root to the  ~/.bagentsetup
# 2.b  Export envirnment variables needed OCDB_PATH, AGENTINPUTTYPE, AGENTSE in the .bagentsetup script
# 2.c  Get a valid token: alien-token-init
# 2.d  Copy certificates to default place -  ~/.agentauth/
#      cp /tmp/*$UID* ~/.agentauth/



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
# test if ~/.bagentsetup exists
if [ ! -f ~/.bagentsetup ]; then
  echo "ERROR: '~/.bagentsetup' init script not found\!"
  exit 1
fi

source ~/.bagentsetup
echo xxxxxxxxxxxxxxxxxxxxxxxxx
echo ENVIRONMENT
echo ALIROOT   =`which aliroot` 
echo    ROOT   =`which root` 
echo xxxxxxxxxxxxxxxxxxxxxxxxx


#chek env
if [ "x$AGENTINPUTTYPE" == "x" ]; then
  echo "WARNING: environment variable 'AGENTINPUTTYPE' falling back to default (castor)\!\!"
  export     AGENTINPUTTYPE=2
fi
if [ "x$OCDB_PATH" == "x" ]; then
  echo "WARNING: environment variable 'OCDB_PATH' falling back to default (local:///afs/cern.ch/alice/tpctest/OCDB)\!\!"
  export     OCDB_PATH=local:///afs/cern.ch/alice/tpctest/OCDB
fi
if [ "x$AGENTSE" == "x" ]; then
  echo "WARNUNG: environment variable 'AGENTSE' falling back to default (ALICE::CERN::SE)\!\!"
  export     AGENTSE=ALICE::CERN::SE
fi

#
#
echo OPERATING SYSTEM
uname -a
#export workdir
export jobhome=`pwd`

#cd to workdir and check for 'macros' and 'out' dir
cd $1
if [ ! -d macros ]; then
  echo "ERROR: no 'macros' dir in '$1'. Does it exist?"
  exit 1
fi
if [ ! -d out ]; then
  echo "ERROR: no 'out' dir in '$1'. Does it exist?"
  exit 1
fi

# check if we are able to write to this dir, the alien SE, can readback from the alien SE and delete stuff
touch agend_alien_tst$HOSTNAME.txt || exit 1
date >> agend_alien_tst$HOSTNAME.txt
#copy to alien
alien_cp -dn agend_alien_tst$HOSTNAME.txt alien:${alien_HOME}agend_alien_tst$HOSTNAME.txt@$AGENTSE

#see if it is there
res=`alien_ls ${alien_HOME}agend_alien_tst$HOSTNAME.txt | grep agend_alien_tst`
if [ "x$res" == "x" ]; then
  echo "ERROR: could not write on storage element '$AGENTSE'"
  exit 1
fi

#copy back
alien_cp -n alien:${alien_HOME}agend_alien_tst$HOSTNAME.txt agend_alien_tst_back.txt 
res=`diff agend_alien_tst$HOSTNAME.txt agend_alien_tst_back.txt`
if [ "x$res" != "x" ]; then
  echo "ERROR: problems reading from storage element '$AGENTSE'"
  exit 1
fi

#try to delete from alien
res=`alien_rm -d ${alien_HOME}agend_alien_tst$HOSTNAME.txt`
if [ "x$res" != "x" ]; then
  echo "ERROR: cannot delete from storage element '$AGENTSE'"
  exit 1
fi

#try to delete local files
res=`rm agend_alien_tst*`
if [ "x$res" != "x" ]; then
  echo "ERROR: cannot delete local files"
  exit 1
fi

echo Current Dir
pwd

aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/AliTPCjobs.cxx+

#alien_cp -d job.list alien:${alien_HOME}job.list@$AGENTSE

