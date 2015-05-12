#!/bin/bash

#  source $ALICE_PHYSICS/../src/PWGPP/rawmerge/validation.sh


     cat <<EOF > exceptions.txt
std::bad_alloc
Segmentation violation
Segmentation fault
Bus error
floating point exception
Killed
busy flag cleared
EOF


cat exceptions.txt


#check log files for presence of exceptions
#if theyre are - exit errorCode


