#!/bin/bash 

# --- Variables needed by base stuff ---------------------------------
base=${ANA_SRC}/scan
datadir=/data/alice/data/pbpb/LHC10h/pass2

# --- Our set-up -----------------------------------------------------
runs="137848 138190"
slCuts="fix=0,0.15"
shCuts="sig=.5,1;xi=1,2;prob=1"
dcCuts="prob=1;sig=.5,1;xi=1,2,5"
strCuts="false"

# --- Help -----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS] [-- [TRAIN OPTIONS]]

Options:
	-h,--help		This help
	-d,--datadir DIR	Top level directory of data
	-r,--runs    RUNS       Runs to analyze
	-2,--no-3		Do not do 3-strip sharing
	-n,--no-act		Do not do any processing, just show
	--			Terminate command line handling 

RUNS is a space separated list of runs 
EOF
}

noact=0
dbg=0
# --- Handle command line --------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)	usage ; exit 0 ;;
	-d|--datadir)   datadir=$2 ; shift ;;
	-D|--debug)     dbg=$2 ; shift ;; 
	-r|--runs)      runs="$2" ; shift ;;
	-2|--no-3)      strCuts="false" ;; 
	-n|--no-act)    noact=1 ;;
	-m|--mc)        mc=1 ;;
	--)             shift ; break;;
	*) echo "$0: Unknown argument $1" >/dev/stderr; exit 1;;
    esac
    shift
done

# --- Source basic setup ---------------------------------------------
. ${base}/baseScan.sh 
debug=$dbg

# --- Loop over everyting --------------------------------------------
fullLoop "$runs" "$strCuts" "$slCuts" "$shCuts" "$dcCuts" $@

# EOF
