#!/bin/bash
/home/perthi/HLT/hlt-current/src/SimpleChainConfig1/MakeTaskManagerConfig.py \
-masternode 0 -taskmandir "/home/perthi/HLT/hlt-current" \
-frameworkdir "/home/perthi/HLT/hlt-current" \
-rundir "/home/perthi/HLT/phos/log/" \
-basesocket 15006 -usessh \
-prestartexec ". /home/perthi/setenv.sh" \
-config "/home/perthi/HLT/phos/config.xml"

