#!/bin/sh 

oper=$1 
if test "x$oper" = "x" ; then oper=init ; fi 

runTrain --class=MyTrain --name=myJob  	\
    --root=v5-33-02b			\
    --aliroot=v5-03-24-AN		\
    --alien=V1.1x			\
    --overwrite				\
    --date=now				\
    --verbose=1				\
    --run=178028			\
    --run=178026			\
    --datadir=/alice/data/2012/LHC12b	\
    --pattern=ESDs/pass1/*/		\
    --mode=grid 			\
    --per-run-merge 			\
    --oper=${oper}
