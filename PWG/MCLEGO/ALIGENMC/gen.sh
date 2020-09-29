#! /bin/bash
# Use SPLIT_MAX_INPUT_FILE_NUMBER in order to get the number of events per job
# STEER energy, tune, packages and pt-hard bins via command line arguments
# of AliGenExtExec

GENERATOR=""
ENERGY=13000
TUNE=""
PACKAGES=""
KTHARDMIN=
KTHARDMAX=

for i in "$@"
do
case $i in
    -g=*|--generator=*)
    ENERGY="${i#*=}"
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
    --kthardmin=*)
    KTHARDMIN="${i#*=}"
    shift
    ;;
    --kthardmax=*)
    KTHARDMAX="${i#*=}"
    shift
    ;;
    *)
    # unknown option
    echo "Unknown option $i, cannot parse"
    ;;
esac
done

# check whether at least the generator is defined
if [ "x$GENERATOR" == "x" ]
then
  echo "Generator needs to be defined"
  exit 1
fi

# Print settings
echo "Generator:             $GENERATOR"
echo "Energy:                $ENERGY"
echo "Seed:                  $ALIEN_PROC_ID"
echo "Number of events:      $SPLIT_MAX_INPUT_FILE_NUMBER" 
echo "Tune:                  $TUNE"
echo "Packages:              $PACKAGES"
echo "Min. pt-hard:          $KTHARDMIN"
echo "Max. pt-hard:          $KTHARDMAX"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv aligenmc::v0.0.5-2)

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