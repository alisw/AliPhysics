#!/bin/bash

set -e
set -x

opts="--branches"
redir="2>/dev/stdout | tee"
redir2=

if test "X$1" == "X--batch" ; then
    opts="${opts} --batch"
    redir=">"
    redir2="2>&1"
    shift 
fi 


./test.sh mc --events=20000 ${opts} $@ > mc.log 2>&1 
./test.sh --events=100000 ${opts} $@ > real.log 2>&1 

