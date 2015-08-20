#!/bin/bash

# Author: Jochen Klein <Jochen.Klein@cern.ch>

eval RUNSCRIPT=$1
FILESTOCOPY="${FILESTOCOPY} ${RUNSCRIPT}"

cat <<EOF > ${RUNSCRIPT}
########################
date
hostname
pwd
########################
EOF

chmod u+x ${RUNSCRIPT}

addtorun() {
    echo $* >> ${RUNSCRIPT}
}
