#!/bin/bash

die() {
    echo $1
    exit 1
}

check_file() {
    [[ ! -f $1 ]] && die "$1 does not exist -> validation failed"
}

check_file galice.root
check_file Kinematics.root
check_file gen.log

exit 0
