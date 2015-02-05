#!/bin/bash

export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
cat <<EOF
=========================================
PATH=$PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH
ROOTSYS=$ROOTSYS
`which root`
ALICE_PHYSICS=$ALICE_PHYSICS
`which aliroot`
`ulimit -a`
`free -m`
Content of directory
`ls -al`
=========================================
EOF

which=$1
run=$2
dir=$3
stage=$4
counter=$5

cat <<EOF
Merge $which Stage $stage sub-job $counter 
  Directory: $dir
  Run:       $run
EOF

if test ! -f ${which}.C ; then 
    echo "No ${which}.C found in current directory" \
	> validation_error.message
    exit 1
fi 
ARG="${which}.C($run,\"$dir\",$stage)"
echo "Running aliroot -b -q -x $ARG"
time aliroot -b -q -x $ARG
exitcode=$?

echo "======== $ARG finished with exit code: $exitcode ========"
echo "############## memory after: ##############"
free -m
if [ $exitcode -ne 0 ]; then
    echo "$ARG exited with code $exitcode" > validation_error.message
fi

exit $exitcode
#
# EOF
#
