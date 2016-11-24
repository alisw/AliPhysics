#!/bin/bash -l

# a script checking the outcome of runtest.sh
# can be used to detect most common errors automatically
# returns/exits with 0 if no problem

# a set of simple rules: (which can be extended)

error=0

value=`grep "RunSimulation  R" sim.log`
if [ -z "${value}" ]; then
  echo "check File not finished"
  let error=error+1
fi

# Pythia errors
value=`grep "Error type*PYEXEC" sim.log`
if [ -z "${value}" ]; then
  echo "Pythia Warning"
# let error=error+1
fi

# reconstruction
value=`grep "RunTracking  R" rec.log`
if [ -z "${value}" ]; then
  echo "RunTracking not finished"
  let error=error+1
fi

# check ESD in check.log
value=`grep "check of ESD was successful" check.log`
if [ -z "${value}" ]; then
  echo "ESD not successful"
  let error=error+1
fi

# check efficiency
line=`grep "I-CheckESD: eff" check.log`
if [ -z "${value}" ]; then
  echo "ESD efficiency missing"
  let error=error+1
else
  eff=`awk '/I-CheckESD: eff/{print $4}' check.log | sed 's/(//'`
  intpart=${eff%.*}
  if (( intpart < 90 )); then
    echo "ESD efficiency less than 90"
    let error=error+1
  fi
fi

exit ${error}
