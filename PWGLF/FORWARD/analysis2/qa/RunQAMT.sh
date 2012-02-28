#!/bin/bash
#
# This script runs the Forward QA for the specified production number
#
# The scripts downloads and runs the single run QA in parallel 
#

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
	-R,--also-results          Also get QAresults.root for each run
	-Q,--qa-number  NUMBER     Custom QA id [$qanumber]
	-p,--production IDENTIFIER Production identifier [$prodfull]
	-P,--pass       IDENTIFIER Pass identifier [$passfull]
	-l,--log-file              Log file output [$redir]

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
qanumber=0
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
# --- Get the index for information ----------------------------------
skip=0
jobu="http://alimonitor.cern.ch/prod/jobs.jsp?t="
jobid=
get_job() 
{ 
    mess 1 "Getting the job information" 
    wget -q ${jobu}${jobid} -O job.html
    local lskip
    let lskip=$skip+1
    p=`grep "/catalogue/index.jsp?path" job.html | head -n $lskip | tail -n 1 | sed -e 's,.*/alice/\(data\|sim\)/\([^<]*\)<.*,\2,' | tr '/' ' '` 
    rm -f job.html
    get_parts $p 
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
other=QAresults.root
files=
path=
numf=0
get_filelist()
{
    mess 3 "Getting file list" 
    
    local datd=data
    local esdd=ESDs/
    case x$prodpost in 
	x_*) ;; 
	x) ;; 
	*)  mess 3 "Assuming simulation output"
	    datd=sim 
	    esdd= 
	    ;; 
    esac
    
    local paid=
    if test "x$passfull" != "x" && test $passno -gt 0 ; then 
	paid=pass${passno}
    fi
    local post=${passpost}
    case x$post in 
	x_*) ;; 
	x) ;; 
	*) post="_${post}" ;; 
    esac

    path=/alice/${datd}/${year}/${prodfull}/
    local search="${esdd}${passpre}${paid}${post}"

    if test $qanumber -gt 0 ; then 
	qapost=`printf "QA%02d" $qanumber` 
	search=`append_to_path "$search" $qapost` 
    fi
    
    search=`append_to_path "$search" $file` 

    cat <<EOF
	Path:			$path
	Search:			$search
EOF
    mess 1 "Getting list of files from AliEn - can take minutes - be patient"
    mess 2 "alien_find ${path} ${search}"
    files=`alien_find ${path} ${search} | grep -v "files found" 2>> ${redir}` 
    for i in $files ; do 
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
    root -l -b  <<EOF >> ${redir} 2>&1 
.L $ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/CheckQAFile.C
CheckQAFile("$1");
.q
EOF
    local ret=$? 
    mess 2 "Check of $1 -> $ret"
    rm -f ${scr}.C 
    return $ret
}

# --- Analyse a file -------------------------------------------------
analyse_file()
{
    local dir=`dirname $1` 
    local inp=`basename $1` 
    local out=`echo $inp | sed 's/trending_/tree_/'` 
    local ret=0
    mess 2 -n "Analysing $inp -> $out ... "

    if test -f $dir/$out ; then 
	mess 2 "exits"
	return 0
    fi

    (cd $dir 
	root -l -b  <<EOF > /dev/stderr
.L $ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/RunFileQA.C
RunFileQA("$inp", "$out");
.q
EOF
	ret=$? 
	mess 2 " -> $ret"
	rm -f ${scr}.C 
    ) 2>> $redir
    return $ret
}

# --- Download a single file -----------------------------------------
also_results=1
analyse_run()
{
    local source=$1 ; shift 
    local store=$1 ; shift 
    local r=$1 ; shift 
    local o=${store}/`basename $file .root`_${r}.root 

    mess 2 -n "$source -> $o ... "
    if test -f $o ; then 
	mess 2 "exists" 
	# sleep 1
    else
	mess 2 -n "copying ... " 
	alien_cp alien:${source} file:${o} >> ${redir} 2>&1 
	fix_perm $o 
	mess 2 "done"
    fi
    if test ! -f $o ; then return 1 ; fi 

    if test $also_results -gt 0 ; then 
	local s=`dirname ${source}`/${other}
	local q=${store}/`basename $other .root`_${r}.root

	mess 2 -n "$s -> $q ... "
	if test -f $q ; then 
	    mess 2 "exists" 
	else
	    mess 2 -n "copying ... " 
	    alien_cp alien:${s} file:${q} >> ${redir} 2>&1 
	    fix_perm $q
	    mess 2 "done"
	fi
    fi

	
    check_file ${o} 
    local ret=$? 
    case $ret in 
	0|2) ;; 
	1|3|4|5|6) return 2 ;; 
    esac

    analyse_file ${o}

    return 0
}

# --- Submit run analysis to background ------------------------------
submit_runs()
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

	printf "%3d/%3d: %s\n" $cur $max $r 
	runs[$counter]=$r

	let counter=$counter+1
	
	analyse_run $i $out $r &
	j=`jobs %% | sed -e 's/^[^0-9]*//' -e 's/[^0-9]*$//'` 
	joblist="$joblist $j"
    done
    
    counter=0
    mess 5 "will wait for jobs $joblist"
    for i in $joblist ; do 
	mess 5 "waiting for $i of $joblist"
	wait %$i
	mess 5 "Analysing ${runs[$counter]} returned $?" 
	let counter=$counter+1
    done
}
    
# --- Analyse each run in turn ---------------------------------------
maxjobs=`grep "processor" /proc/cpuinfo | wc -l`
analyse_runs()
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
	    submit_runs $out $start $max $list 
	    list=
	    queued=0 
	fi
    done
    if test $queued -gt 0 ; then 
	mess 1 "Submitting $queued jobs from $start/$max"
	submit_runs $out $start $max $list 
    fi
}

# --- Copy style -----------------------------------------------------
copy_style()
{
    if test ! -f $style ; then return ; fi 
    rm -f style.css 
    cp $style . 
    fix_perm style.css 
}	

# --- Run the final trending -----------------------------------------
make_trend()
{
    local dir=$1 
    local ret=0
    mess 1 -n "Analysing $dir ... "
    (cd $dir 
	rm -f trend_*_*.html 
	rm -f trend_*_*.pdf
	rm -f trend_*_*.root

	root -l -b <<EOF > /dev/stderr
.L $ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/RunFinalQA.C
RunFinalQA(".");
.q
EOF
	mess 1 -n " ... "
	mess 3 -n "root -l -b -q ${scr}.C "
	root -l -b -q ${scr}.C  > /dev/null 2>&1 
	local ret=$? 
	mess 1 " -> $ret"
	rm -f ${scr}.C 

	# do the index file 
	local idx=`ls trend_*_*.html 2> /dev/null` 
	for i in $idx ; do 
	    mess 1 "Making index.html point to $i" 
	    sed -e 's,index.html,../index.html,' \
		-e "s,!--JOBID--,a target='_blank' href='${jobu}${jobid}'>Job</a," \
		< $i > index.html 
	    cp index.html $i
	done
	
	if test ! -f index.html ; then 
	    echo "No index file found" 
	    ret=1
	else 
	    fix_perm index.html 
	    fix_perm . > /dev/null 2>&1 
	fi

	copy_style
    ) 2>> $redir
    return $ret
}

# --- Make index file ------------------------------------------------
make_index()
{
    local dir=$1   ; shift  
    local title=$1 ; shift
    local desc=$1  ; shift 
    mess 1 "Making index in $dir ($title)"
    
    (cd $dir
	local date=`date` 
	
	rm -f index.html 

	cat <<EOF > index.html
<html>
  <head>
    <title>$title</title>
    <link rel='stylesheet' href='style.css'>
  </head>
  <body>
    <h1>$title</h1>
EOF
	if test ! "x$desc" = "x" ; then 
	    echo "$desc" >> index.html
	fi
	echo "      <ul>" >> index.html
	for i in * ; do 
	    if test ! -d $i ; then continue ; fi 
	    echo "      <li><a href='$i'>$i</a></li>" >> index.html
	done
	echo "      </ul>" >> index.html 
	if test "x$desc" = "x" ; then 
	    echo "      <div class='back'><a href='../'>Back</a></div>" \
		>> index.html
	fi
	cat <<EOF >> index.html    
    <div class='change'>Last update: $date</div>
  </body>
</html>
EOF
	copy_style 
	fix_perm index.html
	fix_perm . > /dev/null 2>&1 
    )
}


# --- Pass command line options --------------------------------------
redir=/dev/null
maxf=-1
top=.
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help) usage ; exit 0 ;; 
	-v|--verbose) let verb=$verb+1 ;; 
	-j|--jobid)   jobid=$2 ; shift ;; 
	-m|--max-files) maxf=$2 ; shift ;; 
	-M|--max-jobs)  maxjobs=$2 ; shift ;;
	-t|--top)       top=$2 ; shift ;;
	-R|--also-results) also_results=1 ;; 
	-Q|--qa-number) qanumber=$2 ; shift ;;
	-p|--production) 
	    prodfull=$2 
	    shift 
	    parse_prod 
	    year=20${prodyear}
	    ;; 
	-P|--pass) 
	    passfull=$2
	    shift
	    parse_pass 
	    ;;
	-l|--log-file) redir= ; shift ;; 
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
    if test "x$prodfull" = "x" || test "x$passfull" = "x" ; then 
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
if test ! "x$passno" = "x" ; then 
    store=${store}/pass${passno}
elif test ! "x$prodpost" = "x" ; then 
    proddir=${proddir}${prodpost}
    store=${proddir}/sim
fi
if test ! "x$qanumber" = "x" && test $qanumber -gt 0 ; then 
    store=${store}_QA${qanumber}
fi
mkdir -p ${top}/$store 
fix_perm ${top}/${proddir}
fix_perm ${top}/$store

if test "x$redir" = "x" ; then 
    redir=${top}/$store/qa.log 
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
EOF
# --- Do a search to find our files ----------------------------------
get_filelist

if test $maxf -gt 0 && test $maxf -lt $numf ; then 
    numf=$maxf 
fi

# --- Now get and analyse each run -----------------------------------
analyse_runs ${top}/$store $numf $files

# --- Now analyse all runs -------------------------------------------
style=$ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/style.css 
make_trend ${top}/$store

# --- Make index files -----------------------------------------------
make_index ${top}/${proddir} ${proddir}
make_index ${top} "QA for the FMD" \
    "For more information see <a href='https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FMDQA'>TWiki pages</a>."

#
# EOF
#
