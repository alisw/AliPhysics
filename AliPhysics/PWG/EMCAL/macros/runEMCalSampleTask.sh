#!/bin/bash

if [ -z "$1" ]; then
	DATATYPE="AOD"
else
	DATATYPE="${1}"
fi

if [ -z "$2" ]; then
	PERIOD="LHC11h"
else
	PERIOD="${2}"
fi

if [ -z "$3" ]; then
	FILELIST="fileLists/files_LHC11h_2_AOD145.txt"
else
	FILELIST="${3}"
fi

if [ -z "$4" ]; then
	NEVENTS="5000"
else
	NEVENTS="${4}"
fi

if [ -z "$5" ]; then
	DEBUGL="0"
else
	DEBUGL="${5}"
fi

root -b -q -x $ALICE_PHYSICS/PWG/EMCAL/macros/runEMCalSampleTask.C\(\""${DATATYPE}"\",\""${PERIOD}"\",\""${FILELIST}"\","${NEVENTS}","${DEBUGL}"\)
