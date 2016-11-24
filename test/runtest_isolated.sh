#!/usr/bin/env bash

# author: Sandro Wenzel ( sandro.wenzel@cern.ch )

# script to execute an AliRoot test in a standalone environment
# provided benefits:
# * the original test directories stay clean
# * we can run the tests in valgrind or other instrumentation or check tools
# * we can perform many tests in various configurations and compare the output

# FIXME: this script is a very rudimentary solution; I'd like to see whether we
# can combine it with the runTests (jenkins) script

# what is this script supposed to do:
# a) copy the content of the test directory to the new location
# (after analysis what needs to be changed)
# b) prepend the aliroot commands with valgrind / perf / igprof or whatever

# Expected input:
# 1 (mandatory) The script expects the test to run (e.g. "gun", or "vmctest/ppbench")
# 2 (mandatory) The script expects a directory which the test is to be run in isolation
# 3 (optional) A command prepended to aliroot (e.g., "valgrind", "valgrind --tool=callgrind", ..)

if [ -z ${1} ] || [ -z ${2} ]; then
  echo "Usage: ${0} testdirectory rundirectory [ prependcommand ]"
  echo " -- Example: ${0} test/gun /tmp/gunrun1/ \"valgrind --tool=callgrind\""
  exit
fi

TESTDIR=$1
TARGETDIR=$2
PREPENDCMD=${3}

# FIXME: check that input makes sense

# copy initial directory
if [ -d "${ALICE_ROOT}" ]; then
  cd ${ALICE_ROOT}
  # check that both the initial and final directories exist

  if [ ! -d "${TESTDIR}" ]; then
    echo "test does not exist; aborting"
    exit
  fi

  if [ -d "$TARGETDIR" ]; then
    echo "run dir already exists; aborting for the sake of security"
    exit
  fi

  if [ ! -d "${TARGETDIR}" ]; then
   # create target directories
    echo "creating run directory ${TARGETDIR}"
   # FIXME: check if dir was given relative or absolute
    mkdir -p ${TARGETDIR}
  fi

  rsync --copy-links -Rr ${TESTDIR}/* ${TARGETDIR}

  FINALDIR="${TARGETDIR}/${TESTDIR}"

  cd $FINALDIR

  # prepend aliroot with prependcommand
  eval "sed -ibak 's/aliroot/${PREPENDCMD} aliroot/' runtest.sh"

  echo "everything setup ... you can now launch runtest.sh in directory ${FINALDIR}"
else
  echo "Please setup the ALICE environment first (expect \$ALICE_ROOT)"
fi
