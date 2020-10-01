#! /bin/bash
# Use SPLIT_MAX_INPUT_FILE_NUMBER in order to get the number of events per job
# STEER energy, tune, packages and pt-hard bins via command line arguments
# of AliGenExtExec

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
echo "Number of events:      $SPLIT_MAX_INPUT_FILE_NUMBER" 
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

#source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv --no-refresh printenv $ALIGENMC_VERSION)


# build command
cmd=$(printf "aligenmc -g %s -E %d -N %d -S %d" $GENERATOR $ENERGY $SPLIT_MAX_INPUT_FILE_NUMBER $ALIEN_PROC_ID)
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