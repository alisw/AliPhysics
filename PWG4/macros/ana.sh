#!/bin/bash

# ana.sh
# 
#
# Created by schutz on 09.04.08.
# Copyright 2008 ALICE. All rights reserved.

#Analysis mode, 0 local, 1 localCAF, 2 PROOF, 3 GRID
export MODE=$1

#Local mode parameters
#Data root directory
export INDIR=$2
#Pattern of data directory name
export PATTERN=Run
#Number of files to analyze
export NEVENT=$3

#GRID mode parameters
#name of xml file
export XML=$2

#Do the analysis
alienroot -b  <<EOF
.L ana.C
ana($MODE)
EOF
