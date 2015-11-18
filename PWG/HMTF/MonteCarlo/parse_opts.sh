# Author: Jochen Klein <Jochen.Klein@cern.ch>

options() {
    declare -f |sed -n -e '/^opt_\S\+ ()\s*$/s/^opt_\(\S\+\).*/\1/p'
}

params() {
    declare -f |sed -n -e '/^par_\S\+ ()\s*$/s/^par_\(\S\+\).*/\1/p'
}

describe_options() {
    echo "options:"
    for opt in $(options); do
	echo "  -$opt         $(opt_$opt)"
    done
    echo
    
    echo "parameters:"
    for par in $(params); do
	echo "  -$par <arg>   $(par_$par)"
    done
    echo
    echo "The values listed for the parameters are what you get with the given command line."
}

parse_args() {
    OPTSTRING=""
    for OPT in $(options); do
	OPTSTRING="${OPTSTRING}${OPT}"
    done
    for PAR in $(params); do
	OPTSTRING="${OPTSTRING}${PAR}:"
    done

    while getopts "${OPTSTRING}" OPT; do
	declare -F opt_${OPT} > /dev/null && opt_${OPT} > /dev/null
	declare -F par_${OPT} > /dev/null && par_${OPT} ${OPTARG} > /dev/null
    done
}

help() {
    printf "%-45s" "$*"
}

setvar() {
    echo "| ${!1}"
    if [ $# -le 1 ]; then
	[ -z ${OPTARG} ] && eval "$1=yes" || eval "$1=${OPTARG}" 
    else
	eval "$1=$2"
    fi
}

show_help() {
    echo "Usage: $(basename $0) [options]"
    echo "${HELP}"
    describe_options
}
