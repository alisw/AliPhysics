#! /bin/bash
####################################################################################
# Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 #
# All rights reserved.                                                             #
#                                                                                  #
# Redistribution and use in source and binary forms, with or without               #
# modification, are permitted provided that the following conditions are met:      #
#     * Redistributions of source code must retain the above copyright             #
#       notice, this list of conditions and the following disclaimer.              #
#     * Redistributions in binary form must reproduce the above copyright          #
#       notice, this list of conditions and the following disclaimer in the        #
#       documentation and/or other materials provided with the distribution.       #
#     * Neither the name of the <organization> nor the                             #
#       names of its contributors may be used to endorse or promote products       #
#       derived from this software without specific prior written permission.      #
#                                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           #
# DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      #
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       #
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     #
####################################################################################
# Run script invoking the aligenmc process in order to generate events
# running an external generator which are read in by aliroot as HepMC 
# input. Furthermore configuring the event generator.
# Steering number of events, energy, tune, packages and pt-hard bins via command line 
# arguments of AliGenExtExec

GENERATOR=""
ENERGY=13000
NEVENTS=
TUNE=""
PACKAGES=""
ALIGENMC_VERSION=""
KTHARDMIN=
KTHARDMAX=
echo $@

for i in "$@"
do
case $i in
    -g=*|--generator=*)
    GENERATOR="${i#*=}"
    shift
    ;;
    -n=*|--nevents=*)
    NEVENTS="${i#*=}"
    shift
    ;;
    -e=*|--energy=*)
    ENERGY="${i#*=}"
    shift
    ;;
    -t=*|--tune=*)
    TUNE="${i#*=}"
    shift
    ;;
    -p=*|--packages=*)
    PACKAGES="${i#*=}"
    shift
    ;;
    -a=*|--aligenmc=*)
    ALIGENMC_VERSION="${i#*=}"
    shift
    ;;
    --pthardmin=*)
    KTHARDMIN="${i#*=}"
    shift
    ;;
    --pthardmax=*)
    KTHARDMAX="${i#*=}"
    shift
    ;;
    *)
    # unknown option
    echo "Unknown option $i, cannot parse"
    shift
    ;;
esac
done

# check whether at least the generator is defined
if [ "x$GENERATOR" == "x" ]
then
  echo "Generator needs to be defined"
  exit 1
fi

if [ "x$NEVENTS" == "x" ]
then
  echo "Number of events / job needs to be defined"
  exit 1
fi

# Print settings
echo "aligenmc:              $ALIGENMC_VERSION"
echo "Generator:             $GENERATOR"
echo "Energy:                $ENERGY"
echo "Seed:                  $ALIEN_PROC_ID"
echo "Number of events:      $NEVENTS" 
echo "Tune:                  $TUNE"
echo "Packages:              $PACKAGES"
echo "Min. pt-hard:          $KTHARDMIN"
echo "Max. pt-hard:          $KTHARDMAX"

# prepare environment
# Add option to load custom aligenmc in local mode
if [ "x$ALIGENMC_VERSION" == "x" ]
then
  ALIGENMC_VERSION="aligenmc::v0.0.5-2"
fi

if [ "x$(echo $ALIGENMC_VERSION | grep :: )" != "x" ]; then
  # cvmfs packages
  source /cvmfs/alice.cern.ch/etc/login.sh
  eval $(alienv printenv $ALIGENMC_VERSION)
else
  # dedicated handling for local builds (for local tests):
  # do not source alienv login script
  # do not refresh the environment
  eval $(alienv --no-refresh printenv $ALIGENMC_VERSION)
fi


# build command
cmd=$(printf "aligenmc -g %s -E %d -N %d -S %d" $GENERATOR $ENERGY $NEVENTS $ALIEN_PROC_ID)
if [ "x$PACKAGES" != "x" ]
then
  cmd=$(printf "%s -p %s" "$cmd" $PACKAGES)
fi
if [ "x$TUNE" != "x" ]
then
  cmd=$(printf "%s -t %s" "$cmd" $TUNE)
fi
if [ "x$KTHARDMIN" != "x" ]
then
  cmd=$(printf "%s -k %s" "$cmd" $KTHARDMIN)
fi
if [ "x$KTHARDMAX" != "x" ]
then
  cmd=$(printf "%s -K %s" "$cmd" $KTHARDMAX)
fi
echo "Running command: \"$cmd\""
eval $cmd
rm gen.hepmc
