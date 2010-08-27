#/bin/bash

# N Events
NEVENTS=10

# TriggerName [Belonging to BarrelMultiplicityTrigger Component]
TRIGGER="H-Barrel_pT_Single-V0001.001"
  
# Path to ESD
ESDPATH="/lustre/alice/jthaeder/data/HEAD_2010-07-09/7TeV/pp_Perugia0/014000"

if [ ! -d ./analysis ] ; then
    mkdir analysis
else
    rm ./analysis/*
fi

# -- Create config CDB object 
aliroot -b -l -q ${ALICE_ROOT}/HLT/trigger/macros/makeTriggerConfigurationObject.C'('\"${TRIGGER}\"')' 2>&1 | tee configLog.log

# -- run chain for ESD file
aliroot -b -l -q ${ALICE_ROOT}/HLT/trigger/macros/HLTTriggerTest.C'('${NEVENTS}','\"${TRIGGER}\"','\"${ESDPATH}\"')' 2>&1 | tee ChainLog.log

