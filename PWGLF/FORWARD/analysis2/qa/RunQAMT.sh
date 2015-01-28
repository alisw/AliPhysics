#!/bin/bash
#
# This script runs the Forwardqq QA for the specified production number
#
# The scripts downloads and runs the single run QA in parallel 
#

# --- Some aux files -------------------------------------------------
if test "X$QA_FWD" = "X" ; then 
    QA_FWD=$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/qa
fi 
style=${QA_FWD}/style.css 
favicon=${QA_FWD}/fmd_favicon.png
logo=${QA_FWD}/fmd_logo.png
script=${QA_FWD}/script.js
topmk=${QA_FWD}/makeIndex.sh

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
print() { 
    local col=$1 ; shift 
    local lvl=$1 ; shift 
    local opt=
    if test $lvl -gt $verb; then return ; fi 
    case $1 in 
	-*) opt=$1 ; shift ;; 
    esac
    echo -e ${opt} "\e[${col}m$*\e[0m"
    # echo -e ${opt} "\e[${col}m$*\e[0m"
}
mess()
{
    print 95 $@ 
}

err() 
{
    echo -e "\e[1mError: \e[91m$*\e[0m" > /dev/stderr 
}
warn()
{
    echo -e "\e[1mWarning: \e[93m$*\e[0m" > /dev/stderr 
}
info() 
{
    print 96 0 $@ 
}
ok() 
{
    print 92 0 "OK"
}

# --- Handling of exit -----------------------------------------------
lock=
handle_exit()
{
    if test "x$lock" = "x" ; then return ; fi 
    if test "x$store" != "x" && test "x${top}" != "x" ; then 
	chmod -R g+rwX ${top}/${proddir} >> ${redir} 2>&1
	chmod -R g+rwX ${top}/$store >> ${redir} 2>&1
    fi
    rm -rf $lock 
}
trap handle_exit EXIT

# --- Handling of errors ---------------------------------------------
last=
handle_err()
{
    err "$last"
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
	err "Another QA process is already running:"
	err "Content of ${lock}:"
	cat $lock > /dev/stderr 
	trap - EXIT
	exit 1
    fi 
    local now=`date` 
    cat <<-EOF > $lock
	Process: $$
	User:    $USER
	Start:   $now
	EOF
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
force=0
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
barrel=0

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
from_local=0
type=data
passid=
mc=0
search=
get_filelist()
{
    mess 3 "Getting file list" 

    local datd=data
    local esdd=ESDs/
    if test ${barrel} -gt 0 ; then
	esdd=
    fi
    if test ${barrel} -gt 1 ; then 	
	file=trending_barrel.root
	other=QAresults_barrel.root
    fi
    if test $mc -gt 0 ; then 
	datd=sim 
	type=sim
	esdd= 
    fi 

    local post=${passpost}
    case x$post in 
	x_*) ;; 
	x) ;; 
	*) post="_${post}" ;; 
    esac

    local paid=
    if test "x${passpre}pass${passno}${post}" != "x$passfull" ; then 
	passpre=
	paid=${passfull}
	post=
    elif echo "$passno" | grep -q -E '^[0-9]*[.]?[0-9]*$' ; then 
	if test "x$passfull" != "x" && test $passno -gt 0 ; then 
	    paid=pass${passno}
	fi
    else
	paid=${passfull}
	passpre=
	post=
    fi
    passid=${paid}
    local spass="${passpre}${paid}${post}"
    if test $mc -gt 0 ; then passid="passMC" ; spass= ; fi 

    search=
    if test "x$path" = "x" ; then 
	# Assume official productions 
	path=/alice/${datd}/${year}/${prodfull}/
	search="${esdd}${spass}"
    else
	search="*"
    fi

    if test $qanumber -gt 0 ; then 
	qapost=`printf "QA%02d" $qanumber` 
	search=`append_to_path "$search" $qapost` 
    fi
    
    search=`append_to_path "$search" $file` 

    cat <<EOF
	Path:			$path
	Search:			$search
EOF
    if test $from_local -lt 1 ; then 

	mess 1 "Get list of files from AliEn - can take minutes - be patient"
	mess 2 "alien_find ${path} ${search}"
	files=`alien_find ${path} ${search} | \
	  grep -v "\(files found\|AND THE\)" 2>> ${redir}` 
    else 
	files=`ls ${top}/${store}/*/input.root | \
	  sed 's,${top}/${store}/,,g'`
    fi
    for i in $files ; do 
	let numf=$numf+1
    done 
    mess 1 -n "Total of $numf files ... "
    ret=$numf
    if test $maxf -lt 0 ; then 
	mess 1 "using all" 
    else
	mess 1 "using $maxf first of these"
	ret=$maxf
    fi
    return $ret
}

# --- Change permissions on files ------------------------------------
fix_perm()
{
    local tgt=$1
    local opt= 
    if test -d $tgt ; then opt="-R" ; fi 
    chmod ${opt} g+rwX $tgt >> /dev/null 2>&1 
    chmod ${opt} o+rX $tgt >> /dev/null 2>&1 
}

# --- Check if a file is OK ------------------------------------------
docheck=1
check_file()
{
    if test $docheck -lt 1 ; then return 0; fi 
    root -l -b  <<EOF >> ${redir} 2>&1 
.L ${QA_FWD}/CheckQAFile.C
CheckQAFile("$1","QA");
.q
EOF
    local ret=$? 
    mess 2 "Check of $1 -> $ret"
    # rm -f ${scr}.C 
    return $ret
}

# --- Analyse a file -------------------------------------------------
analyse_file()
{
    local dir=`dirname $1` 
    local inp=`basename $1` 
    local r=$2
    local out=trending.root 
    # `echo $inp | sed 's/trending_/tree_/'` 
    local ret=0
    mess 2 -n "Analysing $inp -> $out ... "

    if test -f $dir/$out ; then 
	if test $force -lt 1 ; then 
	    mess 2 "exits"
	    return 0
	fi
	rm -f $dir/$out
    fi

    
    mess 3 "runQA.sh '$inp' '$type' '$prodyear' '$prodfull' '$passfull' '$r'"
    # mess 3 "runQA.sh '$inp' '$type' '$prodyear' '$prodfull' '$passid' '$r'"
    (cd $dir 
	for i in QABase QAPlotter QARing QAStructs QATrender ; do 
	    rm -f ${i}*
	    ln -s ../${i}* . 
	done 
	${QA_FWD}/runQA.sh \
	    "$inp" "$type" $prodyear "$prodfull" "$passfull" "$r" > runQA.log 2>&1
	ret=$? ) 
    if test ! -f $dir/trending.root ; then ret=1 ; fi
    mess 2 " -> $ret"
    if test $ret -ne 0 ; then 
	err "Failed to analyse $1"
    fi
    return $ret
}

# --- Download a single file -----------------------------------------
also_results=0
analyse_run()
{
    local source=$1 ; shift 
    local store=$1 ; shift 
    local r=`echo $1 | sed 's/^0*//'` ; shift 
    local rr=`printf %09d $r`
    local o=${store}/${rr}/input.root
    mkdir -p ${store}/${rr}

    mess 2 -n "$source ($store) -> $o ... "
    if test -f $o && test $force -lt 2; then 
	mess 2 "exists" 
	# sleep 1
    else
	rm -f ${o}
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
	if test -f $q && test $force -lt 2 ; then 
	    mess 2 "exists" 
	else
	    rm -rf ${q}
	    mess 2 -n "copying ... " 
	    alien_cp alien:${s} file:${q} >> ${redir} 2>&1 
	    fix_perm $q
	    mess 2 "done"
	fi
    fi

	
    check_file ${o} 
    local ret=$? 
    case $ret in 
	0|2) : ;; 
	1|3|4|5|6) return 2 ;; 
    esac

    analyse_file ${o} ${r}

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

	local r
	if test $from_local -lt 1 ; then 
	    local b=`echo $i | sed -e "s,${path}/*,,"` 
	    if test "x$search" != "x" ; then 
		b=`echo $b | sed -s "s,/*${search},,"`
	    fi
	    r=`echo $b | sed -e "s,/.*,," | sed 's/^0*//'` 
	    # local b=`basename $(dirname $i)`
	    # r=`echo $b | sed 's/^0*//'`
	else
	     r=`basename \`dirname $i\` | sed 's/^0*//'`
	fi

	local m=`printf "%3d/%3d: %s\n" $cur $max $r` 
	info "$m"
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
copy_aliroot_file()
{
    local file=$1 
    if test ! -f $file ; then return ; fi 
    base=`basename $file`
    rm -f $base 
    cp $file $base 
    fix_perm $base
}
copy_style()
{
    copy_aliroot_file $style
}	

# --- Run the final trending -----------------------------------------
variance=1
make_trend()
{
    local dir=$1 
    local ret=0
    info -n "Analysing for trend $dir ... "
    (cd $dir 
	mess 1 "hadd trending.root 000*/trending.root"
	rm -f trending.root 
	hadd -k trending.root 000*/trending.root 
	if test $? -eq 0 && test -f trending.root ; then 
 	  ${QA_FWD}/periodQA.sh trending.root 
	  ret=$?
	else 
	  ret=1
        fi
    ) >>${redir} 2>&1
    if test $ret -ne 0 ; then 
	err "Failed to make trending in $dir"
    else
	ok
    fi
    return $ret
}


# --- Help output ----------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS] -p PRODUCTION [-P PASS]

Options:
  -b,--barrel       MODE       Fetch barrel data              [$barrel]
  -d,--directory    DIR        Search custom AliEn directory  [$path]
  -f,--force                   Force re-run analysis          [$force]
  -h,--help	   	       This help 
  -i,--no-index                Do not make index              [$index]
  -l,--log-file                Log file output                [$redir]
  -L,--local                   Local trending_<X>.root files  [$from_local]
  -m,--max-files    NUMBER     Max number of files to get     [$maxf]
  -M,--max-jobs     NUMBER     Max number of consequtive jobs [$maxjobs]
  -p,--production   IDENTIFIER Production identifier          [$prodfull]
  -P,--pass         IDENTIFIER Pass identifier                [$passfull]
  -Q,--qa-number    NUMBER     Custom QA id                   [$qanumber]
  -R,--also-results            Also get QAresults.root/run    [$also_results]
  -t,--top          DIRECTORY  Output directory               [$top]
  -T,--min-max                 Errors=min/max
  -v,--verbose	  	       Increase verbosity             [$verb]
  -V,--variance                Errors=variance                [$variance]

Production identifiers are of the form LHC11h, LHC11h3, or LHC11h_2. 
Pass identifers are of the form pass2, pass1_HLT, or cpass1.  
If barrel mode>0, then do not assume ESD directory.  
If barrel mode>1, then get trending_barrel.root and QAresults_barrel.root
Option -d is for hand-made QA passes. 
If optiond -d is not specified then official QA passes are assumed.
EOF
}


# --- Parse command line options -------------------------------------
redir=/dev/null
maxf=-1
top=.
index=1
while test $# -gt 0 ; do 
    case $1 in 
	-b|--barrel)       barrel=$2          ; shift ;;
	-d|--directory)	   path=$2	      ; shift ;;
	-f|--force)        let force=$force+1 ;; 
	-h|--help)         usage              ; exit 0 ;; 
	-i|--index)        index=1            ;; 
	--no-index)        index=0            ;;
	-l|--log-file)     redir=             ;; 
	-L|--local)        from_local=1	      ;;
	-m|--max-files)    maxf=$2            ; shift ;; 
	-M|--max-jobs)     maxjobs=$2         ; shift ;;
	-p|--production)   prodfull=$2        ; shift ; parse_prod ;;
	-P|--pass)         passfull=$2        ; shift ; parse_pass ;;
	-Q|--qa-number)    qanumber=$2        ; shift ;;
	-R|--also-results) also_results=1     ;; 
	-t|--top)          top=$2             ; shift ;;
	-T|--min-max)      variance=0	      ;; 
	-v|--verbose)      let verb=$verb+1   ;; 
	-V|--variance)     variance=1         ;;
        -C|--no-check)     docheck=0          ;;
	*) echo "$0: Unknown argument: $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift
done
# === Initial setup ==================================================
# --- Check f AliEn token ------------------------------------------
check_token

# --- Check settings -------------------------------------------------
if test "x$prodfull" = "x" ; then
    err "No production specified" 
    exit 1
fi 
if test "x$prodpost" = "x" && test "x$passfull" = "x" ; then 
    err "No pass specified for non-MC production"
    exit 1
else 
    case x$passfull in
	x)
	    warn "No pass specified, assuming MC"
	    passfull="passMC"
	    mc=1 
	    ;;
	xpassMC*) 
	    warn "MC pass specified"
	    mc=1
	    ;;
	*)  : ;;
    esac
fi 

# --- Construct output -----------------------------------------------
proddir=
passdir=
store=
year=20${prodyear}
if test $mc -gt 0 ; then 
    proddir=${prodfull}
    passdir=$passfull
    store=sim
else 
    proddir=LHC${prodyear}${prodletter}
    store=data
    if test "x$passno" = "x" ; then 
	err "No pass number"
	passdir=$passfull
    else 
	if test "x$passpre" != "x" ; then 
	    passdir=${passpre}
	fi
	passdir=${passdir}pass${passno}
	if test "x$passpost" != "x" ; then 
	    passdir=${passdir}_${passpost}
	fi
	if test "x$remainder" != "x" ; then 
	    passdir=${passdir}/${remainder}
	fi
    fi
fi 
store=${store}/${year}/${proddir}/${passdir}
if test ! "x$qanumber" = "x" && test $qanumber -gt 0 ; then 
    store=${store}_QA${qanumber}
fi
mkdir -p ${top}/$store 
fix_perm ${top}/$store

# --- Check for logging ----------------------------------------------
if test "x$redir" = "x" ; then 
    redir=${top}/$store/qa.log 
    rm -f $redir
    fix_perm $redir
fi

check_lock ${top}/$store

# --- Some friendly information --------------------------------------
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
	Force:                  ${force}
	Use variance:           ${variance}
	Use pre-downloaded:	${from_local}
EOF

# --- Do a search to find our files ----------------------------------
get_filelist
nf=$?
if test $nf -le 0 ; then 
    err "No files to process"
    exit 1
fi
if test $maxf -gt 0 && test $maxf -lt $numf ; then 
    numf=$maxf 
fi

# --- Copy scripts to target and compile -----------------------------
for i in QABase.h QAPlotter.C QARing.h QAStructs.h QATrender.C ; do
    cp ${QA_FWD}/$i ${store}/${i}
    rm -f ${store}/`echo $i | tr '.' '_'`.{so,d}
    fix_perm ${store}/${i}
done
mess 1 "Compiling QATrender.C"
(cd $store && root -l -b <<EOF 
gROOT->LoadMacro("QABase.h++g");
gROOT->LoadMacro("QATrender.C++g");
.q
EOF
) 2>> ${redir}
mess 1 "Compiling QAPlotter.C"
(cd $store && root -l -b <<EOF 
gROOT->LoadMacro("QABase.h++g");
gROOT->LoadMacro("QAPlotter.C++g");
.q
EOF
) 2>> ${redir}
(cd ${store} && for i in *.so *.d ; do fix_perm $i ; done)

# --- Now get and analyse each run -----------------------------------
analyse_runs ${top}/$store $numf $files

# --- Now analyse all runs -------------------------------------------
make_trend ${top}/$store

# --- Make index files -----------------------------------------------
if test $index -gt 0 ; then 
    info -n "Making index ... "
    desc="For more see <a href='http://cern.ch/go/6Bwz'>TWiki</a>"
    $topmk --title "QA for the FMD" \
	--description "$desc" \
	--link \
	--max-depth 4 \
	--frame \
	--output index.html 
    # >> ${redir} 2>&1 
    fix_perm index.html
    copy_aliroot_file $script
    ok
fi 
chmod -R g+rwX ${top}/${proddir} >> ${redir} 2>&1


#
# EOF
#
