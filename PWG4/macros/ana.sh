#!/bin/bash

# ana.sh
# 
#
# Created by schutz on 09.04.08.
# Copyright 2008 ALICE. All rights reserved.

export INDIR=$1
export PATTERN=Run
export NEVENT=$2
export MODE=$3
alienroot -b  <<EOF
.L anaGammaAnalysis.C
anaGammaAnalysis($MODE)
EOF
