#!/bin/bash 

# --- Variables needed by base stuff ---------------------------------
base=${ANA_SRC}/scan
datadir=/data/alice/data/pbpb/LHC10h/pass2

# --- Our set-up -----------------------------------------------------
runs="137848 138190"
methods="sig xi mpv prob"
mpvcuts="0.1 0.3 0.5 0.7"
xicuts="1 2 3 5"
probcuts="1 3 5 7"
sfcuts="0.10 fix 0.7 mpv"
strcuts="true false"
mycuts=

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
	-m,--methods METHODS	Override which methods to scan ($methods)
	-c,--cuts    VALUES	Override cuts (for all methods)
	-s,--sharing CUTS	Override sharing cuts ($sfcuts)
	-n,--no-act		Do not do any processing, just show
	--			Terminate command line handling 

RUNS is a space separated list of runs 

METHODS is a space separated string of one or more of 

  mpv    c = X x Delta_p
  xi     c = Delta_p - X x xi
  sig    c = Delta_p - x x (xi + sigma)
  fix    c = X 
  prob   c: P(Delta<c) < X
  fit    Fit range

where X is the cut parameter.

VALUES is a space separated list of X values. 

Sharing CUTS must consist of 4 parts. 
EOF
}

noact=0

# --- Handle command line --------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)	usage ; exit 0 ;;
	-d|--datadir)   datadir=$2 ; shift ;;
	-r|--runs)      runs="$2" ; shift ;;
	-2|--no-3)      strcuts="false" ;; 
	-m|--methods)   methods="$2" ; shift ;; 
	-c|--cuts)      mycuts="$2"; shift ;;
	-s|--sharing)   sfcuts="$2"; shift ;;
	-n|--no-act)    noact=1 ;;
	--)             shift ; break;;
	*) echo "$0: Unknown argument $1" >/dev/stderr; exit 1;;
    esac
    shift
done

# --- Source basic setup ---------------------------------------------
. ${base}/baseScan.sh 

# --- Loop over everyting --------------------------------------------
for r in $runs ; do                      # Loop over runs 
    echo "=== Now processing run $r" 
    for t in ${strcuts} ; do             # Loop over 3-strip merging 
	echo "===  With allow 3-particle merging $t" 
	for m in $methods ; do               # Loop over methods 
	    echo "===   For method $m" 
	    if test "x$mycuts" = "x" ; then 
		# If we should use the defined cuts 
		case $m in 
		    mpv)    cuts="$mpvcuts" ;; 
		    xi|sig) cuts="$xicuts" ;; 
		    prob)   cuts="$probcuts" ;; 
		    *) echo "Unknown cut type: $m" ; exit 1 ;;
		esac
	    else
		# If cuts are overwritten by command line 
		cuts="$mycuts" 
	    fi
	    for c in $cuts ; do              # Loop over cuts 
		echo "===    With cut value $c" 
		runOne $r $t $sfcuts $c $m $@ # 'sfcuts' expand to 4 args
	    done                             # Loop over cuts 
	done                                 # Loop over methods 
    done                                     # Loop over 3-strip merging 
done                                         # Loop over runs 

# EOF
