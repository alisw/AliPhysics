#/bin/bash

Usage() {
    echo "Usage: ${0##*/} <path_to_merged> [<path_to_merged1> ...] "
    exit
}

[[ $# -lt 1 ]] &&  Usage && exit

script="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procRunResid.sh"

[[ ! -f "$script" ]] && echo "did not find script $script" && exit 

for line in "$@"
do
    run=`echo $line | perl -lne 'm|(000)(\d{6})|g&&print "$2"'`
    echo "Processing run $run from path $line"
    $script "$line" "$run" >& proc"$run"".log"
done


