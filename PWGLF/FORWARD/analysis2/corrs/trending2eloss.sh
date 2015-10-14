#!/bin/bash

fwd=$ALICE_PHYSICS/PWGLF/FORWARD/analysis2
if test "x$ANA_SRC" ; then fwd=$ANA_SRC ; fi
noserial=1
force=0
refit=0

# --- Help output ----------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS] -d [DIRECTORY] -p [PATTERN]

Options:
	-d,--directory  DIRECTORY  Directory to seach [$dir]
	-f,--force      MASK       Force download and reprocess 
	-h,--help		   This help 
	-m,--max-files  NUMBER     Max number of files to get [$maxf]
	-M,--max-jobs   NUMBER     Max number of consequtive jobs [$maxjobs]
	-n,--no-action             Do a dry run
	-p,--pattern    FILENAME   File name pattern [$pattern]
	-R,--refit                 Run full re-fitting [$refit]
	-t,--top        DIRECTORY  Output directory [$top]
	-v,--verbose		   Increase verbosity [$verb]

Force mask is a bit pattern of 

  0x1: Force download 
  0x2: Force extract 
  0x4: Force draw 
EOF
}

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
# --- Extract from trending.root -------------------------------------
extract()
{
    local outd=$1 ; shift 
    local outf=$1 ; shift    
    local corr="${outd}/fmd_corrections.root"
    if test -f ${outd}/${outf} ; then
	# --- Check if corr exists, or is newer than downloaded file -
	if test ! -f ${corr} || test ${outd}/${outf} -nt ${corr} ; then
	    rm -f ${corr}
	    (cd ${outd} && \
		    root -l -b -q ${fwd}/corrs/Trending2ELoss.C+\(\"${outf}\"\))
	fi
    else
	echo "${outd}/${outf} does not exists"
	return;
    fi

    # --- If the extraction failed, go on ----------------------------
    if test -f ${outd}/bad.root ; then
	return
    fi 
    plot=${outd}/forward_elossfits.pdf
    if test -f ${outd}/fmd_corrections.root ; then
	# --- Check if corr exists, or is newer than downloaded file -
	if test ! -f ${plot} || test ${corr} -nt ${plot} ; then
	    rm -f ${plot}
	    (cd $outd && \
		    root -l -b -q $fwd/corrs/DrawCorrELoss.C\(0,\"$outf\"\)) \
		> /dev/null 2>&1 
	fi
    else
	echo "$corr does not exists"
    fi 
}

# --- Full refit -----------------------------------------------------
full_refit()
{
    local outd=$1 ; shift 
    local outf=$1 ; shift
    local reft="forward_eloss.root"
    local corr="${top}/${outd}/fmd_corrections.root"

    # --- Check if we should refit -----------------------------------
    if test -f ${outd}/${outf} ; then
	mess 5 "Got ${outd}/${outf}"
	# --- Check if corr exists, or is newer than downloaded file -
	if test ! -f ${outd}/${reft} || \
		test ${outd}/${outf} -nt ${outd}/${reft} ; then
	    rm -f ${outd}/${reft}
	    mess 1 "Refitting ${outd}/${outf} -> ${outd}/${reft}"
	    (cd ${outd} && \
		    root -l -b -q \
			 ${fwd}/corrs/RerunELossFits.C\(0,\"${outf}\",1,\"${reft}\",0x2\) > refit.log 2>&1)
	fi
    else
	echo "${outd}/${outf} does not exists"
	return;
    fi

    # --- If the refit failed, bail out ------------------------------
    if test ! -f ${outd}/${reft} ; then
	mess 3 "${outd}/${reft} does not exists"
	return;
    fi

    # --- Now try to extract -----------------------------------------
    (cd $outd && root -l -b -q $fwd/corrs/ExtractELoss.C >> refit.log 2>&1)
    
    # --- If the extraction failed, go on ----------------------------
    if test ! -f ${corr} || test -f ${outd}/bad.root ; then
	mess 1 "Missing correction or bad"
	rm -rf ${corr}
	return
    fi 
    plot=${outd}/forward_eloss.pdf
    if test -f ${corr} ; then
	# --- Check if corr exists, or is newer than downloaded file -
	if test ! -f ${plot} || test ${corr} -nt ${plot} ; then
	    rm -f ${plot}
	    (cd $outd && \
		    root -l -b -q $fwd/corrs/DrawCorrELoss.C\(0,\"${reft}\"\))\
		> /dev/null 2>&1 
	fi
    else
	echo "$corr does not exists"
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
    mess 10 "_download_file src=$1 out=$2 number=$3 current=$4 " \
	 "max=$5 redir=$redir noact=$6 friends=$friends" 
    outd=`printf %09d $3`
    outf=`basename $pattern`
    _download_file "$1" "$2" "$3" "$4" "$5" "$redir" "$6" "1"
    mess 10 "Will process $outd/${outf} ($1)"

    if test $refit -gt 0 ; then
	full_refit "${top}/${outd}" "${outf}"
    else 
	extract "${top}/${outd}" "${outf}"
    fi 
    
}

# --- Build our code -------------------------------------------------
buildCode()
{
    echo "Building code"
    root.exe -l -b <<EOF
gSystem->AddIncludePath("-I${fwd}");
gROOT->Macro("$fwd/scripts/LoadLibs.C");
gROOT->LoadMacro("${fwd}/corrs/Trending2ELoss.C+g");
Info("", "Build the code");
.q
EOF
}

# --- Source library -------------------------------------------------
. ${fwd}/scripts/libGridDownload.sh

maxjobs=`grep "processor" /proc/cpuinfo | wc -l`
# --- Pass command line options --------------------------------------
redir=/dev/null
maxf=-1
top=.
dir=
pattern=trending.root
friends=0
period=
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)       usage            ; exit 0 ;; 
	-v|--verbose)    let verb=$verb+1 ;; 
	-m|--max-files)  maxf=$2          ; shift ;; 
	-M|--max-jobs)   maxjobs=$2       ; shift ;;
	-t|--top)        top=$2           ; shift ;;
	-l|--log-file)   redir=           ;; 
	-d|--directory)  dir=$2           ; shift ;;
	-n|--no-action)  noact=1          ;;
	-p|--pattern)    pattern=$2       ; shift ;;
	-f|--force)      force=$2         ; shift ;;
	-R|--refit)      refit=1          ; shift ;;
	-L|--period)     period=$2        ; shift ;;
	-P|--pass)       pass=$2          ; shift ;;
	*) echo "$0: Unknown argument: $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift
done
if test "x$dir" = "x" && test "x$period" != "x" ; then
    year=`echo $period | sed 's/LHC\(..\).*/\1/'`
    dir="/alice/data/20${year}/${period}"
    pattern="/${pass}"
    case $pass in
	cpass*_pass*) pattern="${pattern}/trending_barrel.root" ;;
	pass*)        pattern="${pattern}/trending.root" ;;
	*)            pattern="${pattern}/trending.root" ;;
    esac
    top=${period}_${pass}
fi 


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

# --- Build our code -------------------------------------------------
buildCode

# --- Friendly message -----------------------------------------------
cat <<EOF
	Directory:		$dir 
	Path:                   $path
	Pattern:		$pattern 
	Output directory:	${top}
	Lock file:		${lock}
	Log:                    ${redir}
	Max # Files:            ${maxf}
	Max # consequtive jobs: ${maxjobs}
	Force flags:            ${force}
	Full refit:             ${refit}
EOF

# --- In case of forces, we remove the files here --------------------
if test $((${force} & 0x1)) -ne 0 ; then
    read -p "Will remove downloaded trending files, are you sure?" -n 1 ans
    case x$ans in
	xy|xY)
	    find ${top} -name "trending.root" | xargs echo
	    force=$(($force | 0x2))
	    ;;
	*)     echo "Removing nothing" ;;
    esac
fi
if test $((${force} & 0x2)) -ne 0 ; then
    mess 1 "Will remove extracted corrections files"
    find ${top} \
	 -name "fmd_corrections.root" -or \
	 -name "forward_eloss.root" -or \
	 -name "bad.root" -or \
	 -name "diagnostics.root" -or \
	 -name "Upload.C" \
	 -name "refit.log" \
	| xargs rm -f
    force=$(($force | 0x4))
fi
if test $((${force} & 0x4)) -ne 0 ; then
    mess 1 "Will remove summary plots"
    find ${top} -name "forward_eloss.pdf" -or -name "forward_elossfits.pdf" \
	| xargs rm -f
fi


# --- Do a search to find our files ----------------------------------
_get_file_list "$dir" "$pattern" $maxf

mess 3 "Number of files: '$numf'"

if test $maxf -ge 0 && test $maxf -lt $numf ; then 
    numf=$maxf 
fi

# --- Now get and analyse each run -----------------------------------
download_files ${top}/ $numf $maxjobs $files

# --- Merge all made corrections into a single file ------------------
if test ! -d ${top} ; then exit 0 ; fi 
corrs=`find ${top} -name "fmd_corrections.root"`  
for c in ${corrs} ; do
    mess 1 "Merging ${c}"
    root -l -b -q ${fwd}/corrs/Upload.C\(\"${top}\",\"${c}\"\) \
	 >> upload.log 2>&1 
done

    
#
# EOF
#
