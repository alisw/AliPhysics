#!/bin/bash
#
# This script runs the Forward QA for the specified production number
#
# The scripts downloads and runs the single run QA in parallel 
#
noact=0

# --- Help output ----------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS] -d [DIRECTORY] -p [PATTERN]

Options:
	-d,--directory  DIRECTORY  Directory to seach [$dir]
	-f,--friends               Also download friends 
	-h,--help		   This help 
	-j,--jobid      JOBID      Download from a particular job
	-l,--log-file              Log file output [$redir]
	-m,--max-files  NUMBER     Max number of files to get [$maxf]
	-M,--max-jobs   NUMBER     Max number of consequtive jobs [$maxjobs]
	-n,--no-action             Do a dry run
	-p,--pattern    FILENAME   File name pattern [$file]
	-t,--top        DIRECTORY  Output directory [$top]
	-v,--verbose		   Increase verbosity [$verb]
EOF
}


# --- Source library -------------------------------------------------
. $ALICE_PHYSICS_SOURCE/PWGLF/FORWARD/analysis2/scripts/libGridDownload.sh

# --- Download the files ---------------------------------------------
download_files()
{
    local out=$1 ; shift 
    local max=$1 ; shift
    local maxjobs=$1 ; shift 
    
    local queued=0
    local counter=0
    local start=0
    local list=
    while test $# -gt 0 ; do 
	if test $counter -ge $max ; then 
	    break;
	fi

	list="$list $1" 
	shift 
	let queued=$queued+1 
	let counter=$counter+1 

	if test $queued -eq 1 ; then 
	    start=$counter
	fi
	
	if test $queued -ge $maxjobs ; then 
	    mess 1 "Submitting $queued jobs from $start/$max"
	    # $1: Output directory
	    # $2: Sarting off-set
	    # $3: Maximum jobs
	    # $4: Maximum number of files
	    # $5: noact
	    mess 2 "Submit_jobs out=$out start=$start maxjobs=$maxjobs " \
		 "max=$max noact=$noact list=$list" 
	    submit_jobs "$out" "$start" "$maxjobs" "$max" "$noact" $list 
	    list=
	    queued=0 
	fi
    done
    if test $queued -gt 0 ; then 
	mess 1 "Submitting $queued jobs from $start/$max"
	submit_jobs $out $start $maxjobs $max $noact $list 
    fi
}
# --- Download a single file -----------------------------------------
#  $1: Source file
#  $2: Output directory
#  $3: inherit number 
#  $4: Current count
#  $5: Maximum
#  $6: noact

download_file()
{
    #  $1: Source file
    #  $2: Output directory
    #  $3: Number
    #  $4: Current count
    #  $5: Maximum
    #  $6: redirect
    #  $7: noact
    #  $8: multi 
    mess 2 "_download_file src=$1 out=$2 number=$3 current=$4 " \
	 "max=$5 redir=$redir noact=$6 friends=$friends" 
    _download_file "$1" "$2" "$3" "$4" "$5" "$redir" "$6" "$friends"
    if test $friends -gt 0 ; then
	case "x$1" in
	    x*AliESDs*) tmp=`echo $1 | sed 's/AliESDs/AliESDfriends/'` ;;
	    x*) tmp= ;;
	esac
	if test "X$tmp" != "X" ; then
	    _download_file "$tmp" "$2" "$3" "$4" "$5" "$redir" "$6" "$friends"
	fi
    fi
}


maxjobs=`grep "processor" /proc/cpuinfo | wc -l`
# --- Pass command line options --------------------------------------
redir=/dev/null
maxf=-1
top=.
dir=
pattern=
friends=0
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)       usage            ; exit 0 ;; 
	-v|--verbose)    let verb=$verb+1 ;; 
	-j|--jobid)      jobid=$2         ; shift ;; 
	-m|--max-files)  maxf=$2          ; shift ;; 
	-M|--max-jobs)   maxjobs=$2       ; shift ;;
	-t|--top)        top=$2           ; shift ;;
	-l|--log-file)   redir=           ;; 
	-d|--directory)  dir=$2           ; shift ;;
	-n|--no-action)  noact=1          ;;
	-p|--pattern)    pattern=$2       ; shift ;;
	-f|--friends)    friends=1        ;;
	*) echo "$0: Unknown argument: $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift
done
# --- Initial setup --------------------------------------------------
# First, check we have a valid AliEn token, then retrieve the job 
# information to parse out the location of the files we need, and 
# finally make our output directory and check the lock 
check_token

if test "x$dir" = "x" ; then
    echo "No directory specified" > /dev/stderr
    exit 1
fi
if test "x$pattern" = "x" ; then
    echo "No pattern specified" > /dev/stderr
    exit 1
fi
path=$dir

check_lock ${top}/

cat <<EOF
	Directory:		$dir 
	Path:                   $path
	Pattern:		$pattern 
	Output directory:	${top}
	Lock file:		${lock}
	Log:                    ${redir}
	Max # Files:            ${maxf}
	Max # consequtive jobs: ${maxjobs}
EOF
# --- Do a search to find our files ----------------------------------
_get_file_list "$dir" "$pattern" $maxf

mess 3 "Number of files: '$numf'"

if test $maxf -ge 0 && test $maxf -lt $numf ; then 
    numf=$maxf 
fi

# --- Now get and analyse each run -----------------------------------
download_files ${top}/ $numf $maxjobs $files

#
# EOF
#
