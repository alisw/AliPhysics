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
Usage: $0 [OPTIONS] -j [JOBID]
       $0 [OPTIONS] -p PRODUCTION -P PASS 

Options:
	-h,--help		   This help 
	-j,--jobid	NUMBER	   The master job id of the production [$jobid]
	-v,--verbose		   Increase verbosity [$verb]
	-m,--max-files  NUMBER     Max number of files to get [$maxf]
	-M,--max-jobs   NUMBER     Max number of consequtive jobs [$maxjobs]
	-t,--top        DIRECTORY  Output directory [$top]
	-p,--production IDENTIFIER Production identifier [$prodfull]
	-P,--pass       IDENTIFIER Pass identifier [$passfull]
	-l,--log-file              Log file output [$redir]
	-f,--file       FILENAME   Name of files to get [$file]
	-r,--run        NUMBER     Run number [$run]

Note the option -j and the options -p and -P are mutually exclusive,
The option -Q is only used if the options -p and -P are given.
Production identifiers are of the form LHC11h, LHC11h3, or LHC11h_2. 
Pass identifers are of the form pass2, pass1_HLT, or cpass1.
EOF
}

# --- Check AliEn token ----------------------------------------------
check_token()
{
    uid=`id -u`
    genv_file=/tmp/gclient_env_${uid}
    
    if test ! -f ${genv_file} ; then 
	echo "No such file: ${genv_file}, please do alien-token-init" \
	    >/dev/stderr
	exit 1
    fi
    . ${genv_file}
    alien-token-info | grep -q "Token is still valid"
    if test $? -ne 0 ; then 
	echo "Token not valid, please re-new" > /dev/stderr 
	exit 1
    fi
}

# --- Diagnostics output ---------------------------------------------
verb=0
mess()
{
    if test $1 -gt $verb ; then return ; fi 
    shift
    echo $*
}

# --- Handling of exit -----------------------------------------------
lock=
handle_exit()
{
    if test "x$lock" = "x" ; then return ; fi 
    rm -rf $lock 
}
trap handle_exit EXIT

# --- Handling of errors ---------------------------------------------
last=
handle_err()
{
    echo "Error: $last" 
    exit 1
}
enable_trap()
{
    trap handle_err ERR
}
disable_trap()
{
    trap - ERR
}

# --- Check the lock -------------------------------------------------
check_lock()
{
    local loc=$1 
    lock=$loc/.lock

    if test -f $lock ; then 
	echo "Another QA process is already running:" > /dev/stderr 
	echo "Content of ${lock}:" > /dev/stderr 
	cat $lock > /dev/stderr 
	trap - EXIT
    exit 1
    local now=`date` 
    cat <<EOF > $lock
Process: $$
User:    $USER
Start:   $now
EOF
    fi
}

# --- Parse production information -----------------------------------
parse_prod()
{
    prodyear=`echo $prodfull | sed 's/LHC\(..\).*/\1/'` 
    prodletter=`echo $prodfull | sed "s/LHC${prodyear}\(.\).*/\1/"` 
    prodpost=`echo $prodfull | sed "s/LHC${prodyear}${prodletter}//"` 
}

parse_pass()
{
    passno=`echo $passfull | sed 's/.*pass\([0-9]*\).*/\1/'`  
    passpost=`echo $passfull | sed -n "s/.*pass${passno}_//p"` 
    passpre=`echo $passfull | sed -n "s/pass.*//p"` 
}
# --- Extract parts from the found path ------------------------------
year=0
passfull=
passno=0
passpost=
passpre=
prodfull=
prodyear=0
prodletter=
prodpost=
remainder=
get_parts()
{
    mess 1 "Parsing information from job $@" 
    year=$1 ; shift 
    prodfull=$1 ; shift 
    local lrunn=$1 ; shift
    local ltype=$1 ; shift 
    passfull=$1 ; shift 
    remainder=$1 

    mess 10 "year=$year" 
    mess 10 "prodfull=$prodfull" 
    mess 10 "lrunn=$lrunn" 
    mess 10 "ltype=$ltype" 
    mess 10 "passfull=$passfull" 
    mess 10 "remainder=$remainder"

    if test "x$passfull" = "x" ; then 
	remainder=
	passfull=$ltype 
    fi
    case x$passfull in 
	*pass*) ;; 
	*) remainder=$passfull 
	    passfull= 
	    ;;
    esac
    parse_pass 
    parse_prod

    case x$remainder in 
	xQA*) qanumber=`echo $remainder | sed 's/QA//'`  ;; 
	*) ;; 
    esac

}
# --- Append path element --------------------------------------------
append_to_path()
{
    local tmp=$1 ; shift 
    local add=$1
    if test "x$tmp" != "x" ; then tmp="${tmp}/" ; fi 
    echo ${tmp}${add}
}

# --- Get a list of files to get -------------------------------------
file=trending.root
files=
path=
numf=0
mc=0
get_filelist()
{
    mess 3 "Getting file list" 
    
    local datd=data/
    local esdd=ESDs/    
    local yerd=$year/
    case x$prodpost in 
	x_*) ;; 
	x) ;; 
	*)  mess 3 "Assuming simulation output"
	    datd=sim/ 
	    esdd= 
	    if test $year -lt 2013 ; then 
		yerd=
	    fi
	    mc=1
	    passfull=
	    passpost=
	    passno=0
	    passpre=
	    ;; 
    esac
    
    local paid=
    if test "x$passfull" != "x" && \
	test $passno -gt 0 && \
	test $mc -lt 1; then 
	paid=pass${passno}
    fi
    local post=
    if test $mc -lt 1; then 
	post=${passpost}
	case x$post in 
	    x_*) ;; 
	    x) ;; 
	    *) post="_${post}" ;; 
	esac
    fi

    path=/alice/${datd}${yerd}/${prodfull}/
    local search="$file"

    if test "x$run" != "x" ; then 
	case x$datd in 
	    xsim/) rund=$run ;; 
	    *)     rund=`printf %09d $run` ;; 
	esac
	path=$path/$rund/
    fi
    path=${path}${esdd}${passpre}${paid}${post}

    # search=`append_to_path "$search" $file` 

    cat <<EOF
	Path:			$path
	Search:			$search
EOF
    mess 1 "Getting list of files from AliEn - can take minutes - be patient"
    mess 2 "alien_find ${path} ${search}"
    files=`alien_find ${path} ${search} | grep -v "\(does not\|files found\)" 2>> ${redir}` 
    rm -f .list
    for i in $files ; do 
	echo $i >> .list 
	let numf=$numf+1
    done 
    mess 1 -n "Total of $numf files ... "
    if test $maxf -lt 0 ; then 
	mess 1 "using all" 
    else
	mess 1 "using $maxf first of these"
    fi
}

# --- Change permissions on files ------------------------------------
fix_perm()
{
    if test ! -f $1 ; then return ; fi 
    chmod g+rwX $1
    chmod o+rX $1
}
# --- Check if a file is OK ------------------------------------------
docheck=1
check_file()
{
    
    if test $docheck -lt 1 ; then return 0; fi 
    case $1 in 
	*.root) ;;
	*.zip) return 0 ;;
    esac
	    
    root -l -b  <<EOF >> ${redir} 2>&1 
.L $ALICE_PHYSICS/PWGLF/FORWARD/analysis2/qa/CheckQAFile.C
CheckQAFile("$1");
.q
EOF
    local ret=$? 
    mess 2 "Check of $1 -> $ret"
    rm -f ${scr}.C 
    return $ret
}

# --- Download a single file -----------------------------------------
download_file()
{
    local source=$1 ; shift 
    local store=$1 ; shift 
    local r=$1 ; shift 
    local o=${store}/ 
    case $file in 
	*.root)
	    o=${o}`basename $file .root`_`printf %04d ${r}`.root 
	    ;;
	*.zip)
	    local d=`printf %04d ${r}` 
	    mkdir -p ${o}/${d}
	    o=${o}/${d}/${file}
	    ;;
    esac
    printf "%4d/%4d: %20s -> %20s ..." $cur $max $source $o 

    mess 2 -n "$source -> $o ... "
    if test -f $o ; then 
	printf "exists\n"
	# mess 2 "exists" 
	# sleep 1
    else
	printf "copying\n"
	# mess 2 -n "copying ... " 
	if test $noact -lt 1 ; then 
	    alien_cp alien:${source} file:${o} >> ${redir} 2>&1 
	    fix_perm $o 
	else 
	    sleep 1
	fi
	mess 2 "done"
    fi
    if test $noact -gt 0 ; then return 0 ; fi 
    if test ! -f $o ; then return 1 ; fi 
	

    case $o in 
	*.root) 
	    check_file ${o} 
	    local ret=$? 
	    case $ret in 
		0|2) ;; 
		1|3|4|5|6) return 2 ;; 
	    esac
	    ;;
	*.zip)
	    d=`dirname $o` ;
	    b=`basename $o` ; 
	    mess 3 "Unzipping $b in $d"
	    if ! unzip -n -qq -l $o > /dev/null 2>&1 ; then
		mess 1 "Bad zip file: $o" 
		rm -rf $d 
	    fi
	    # (cd $d && unzip -n -qq $b)
	    ;;
    esac
    # analyse_file ${o}

    return 0
}

# --- Submit run analysis to background ------------------------------
submit_jobs()
{
    local out=$1 ; shift
    local sta=$1 ; shift 
    local max=$1 ; shift
    
    local joblist=
    local counter=0
    mess 5 "Submitting $maxjobs jobs from $sta/$maxf" 
    for i in $@ ; do 
	let cur=$sta+$counter

	local b=`echo $i | sed -e "s,${path},,"` 
	local r=`echo $b | sed -e "s,/.*,,"` 


	let counter=$counter+1

	download_file $i $out $cur &
	j=`jobs %% | sed -e 's/^[^0-9]*//' -e 's/[^0-9]*$//'` 
	joblist="$joblist $j"
    done
    
    counter=0
    mess 5 "will wait for jobs $joblist"
    for i in $joblist ; do 
	mess 5 "waiting for $i of $joblist"
	wait %$i
	let counter=$counter+1
    done
}
    
# --- Analyse each run in turn ---------------------------------------
maxjobs=`grep "processor" /proc/cpuinfo | wc -l`
download_files()
{
    local out=$1 ; shift 
    local max=$1 ; shift 
    
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
	    submit_jobs $out $start $max $list 
	    list=
	    queued=0 
	fi
    done
    if test $queued -gt 0 ; then 
	mess 1 "Submitting $queued jobs from $start/$max"
	submit_jobs $out $start $max $list 
    fi
}

maxjobs=`grep "processor" /proc/cpuinfo | wc -l`
# --- Pass command line options --------------------------------------
redir=/dev/null
maxf=-1
top=.
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)       usage            ; exit 0 ;; 
	-v|--verbose)    let verb=$verb+1 ;; 
	-j|--jobid)      jobid=$2         ; shift ;; 
	-m|--max-files)  maxf=$2          ; shift ;; 
	-M|--max-jobs)   maxjobs=$2       ; shift ;;
	-t|--top)        top=$2           ; shift ;;
	-l|--log-file)   redir=           ;; 
	-f|--file)       file=$2          ; shift ;;
	-n|--no-action)  noact=1          ;;
	-r|--run)        run=$2           ; shift ;;
	-p|--production) prodfull=$2      ; shift ; parse_prod 
	    year=20${prodyear}            ;; 
	-P|--pass)       passfull=$2      ; shift ; parse_pass    ;;
	*) echo "$0: Unknown argument: $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift
done
# --- Initial setup --------------------------------------------------
# First, check we have a valid AliEn token, then retrieve the job 
# information to parse out the location of the files we need, and 
# finally make our output directory and check the lock 
check_token

if test ! "x$jobid" = x ; then 
    if test ! "x$prodfull" = "x" || test ! "x$passfull" = "x" ; then 
	cat <<EOF > /dev/stderr
Option -j ${jobid} and options -p and -P are mutually exclusive 
EOF
	exit 1
    fi
    get_job
else 
    if test "x$prodfull" = "x" && test "x$passfull" = "x" ; then 
	cat<<EOF > /dev/stderr
When specifying prodcution and/or pass both options -p and -P _must_ 
be specified. 
EOF
	exit 1
    elif test ! "x$jobid"  = "x" ; then 
	cat <<EOF > /dev/stderr
Option -j and options -p ${prodfull} and -P ${passfull} are mutually exclusive 
EOF
	exit 1
    fi
fi	

proddir=LHC${prodyear}${prodletter}
store=${proddir}
if test ! "x$prodpost" = "x" ; then   
    proddir=${proddir}${prodpost}
    store=sim/${proddir}
elif test "x$passfull" != "x" && test $passno -gt 0 ; then 
    store=${store}/pass${passno}
fi
if test "x$run" ; then 
    store=${store}/`printf %06d $run`
fi
mkdir -p ${top}/$store 
fix_perm ${top}/${proddir}
fix_perm ${top}/$store

if test "x$redir" = "x" ; then 
    redir=${top}/$store/download.log 
    rm -f $redir
    fix_perm $redir
fi

check_lock ${top}/$store

cat <<EOF
	Year:			$year
	Production:		$prodfull 
	  Year:			$prodyear
	  Letter:		$prodletter
	  Suffix:		$prodpost
	Pass:			$passfull
	  Number:		$passno
	  Prefix:		$passpre
	  Postfix:		$passpost
	Remainder		$remainder
	  QA number		$qanumber
	Output directory:	${store}
	Lock file:		${lock}
	Log:                    ${redir}
	Run:                    ${run}
EOF
# --- Do a search to find our files ----------------------------------
get_filelist

if test $maxf -ge 0 && test $maxf -lt $numf ; then 
    numf=$maxf 
fi

# --- Now get and analyse each run -----------------------------------
download_files ${top}/$store $numf $files

#
# EOF
#
