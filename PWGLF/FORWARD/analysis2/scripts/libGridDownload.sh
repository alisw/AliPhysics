#!/bin/bash
#
# Library of some download functions 

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
noact=0
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
    local nam=$2 

    if test -f $lock ; then 
	echo "Another $nam process is already running:" > /dev/stderr 
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

# --- Append path element --------------------------------------------
append_to_path()
{
    local tmp=$1 ; shift 
    local add=$1
    if test "x$tmp" != "x" ; then tmp="${tmp}/" ; fi 
    echo ${tmp}${add}
}

# --- Get list of files ----------------------------------------------
# Output is stored in .list
spath=
_get_file_list()
{
    local path=$1
    local search=$2
    local maxf=$3
    if test x$maxf = x ; then maxf=-1 ; fi
    spath=$path
    
    mess 1 "Getting list of files from AliEn - can take minutes - be patient"
    mess 2 "alien_find ${path} ${search}"
    files=`alien_find ${path} ${search} | grep -v "\(does not\|files found\)" 2>> ${redir}` 
    rm -f .list
    numf=0
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
# --- Make output file name ------------------------------------------
#
# $1: Source file 
# $2: Top of store
# $3: run
# $4: current file number
# $5: Maximum number of files
# $6: Multi-file download  
_format_out_name()
{
    local source=$1 ; shift 
    local store=$1  ; shift 
    local r=$1      ; shift
    local cur=$1    ; shift
    local max=$1    ; shift
    local multi=$1  ; shift 
    local o=`echo ${store} | sed 's/\/*$/\//'`
    local file=`basename $source`
    mess 10 "Format name from source=$source store=$store r=$r " \
	 "cur=$cur max=$max multi=$multi" > /dev/stderr 
    case x$r in
	x*.*)
	    rr=`echo $r | sed -e 's/^..0*//' -e 's/...\..*//'`
	    mm=`echo $r | sed -e 's/.*\(...\)\..*/\1/' -e 's/^0//'`
	    lm=`echo $r | sed -e 's/.*\.//' -e 's/^0//'`
	    r=$rr
	    cur=$mm
	    max=$lm
	    # cur=$((${mm}*10000+${lm}))
	    ;;
	x)  r=0
	    ;;
    esac
    if test "x$noserial" = "x" || test $noserial -gt 0 ; then
	sub=`printf %09d ${r}`
    else 
	sub=`printf %09d_%04d_%04d ${r} ${cur} ${max}`
    fi
    tmp=$file 
    if test $multi -le 0 ; then 
	case $file in 
	    *.root)
		tmp=`basename $file .root`_${sub}.root
		sub=		
		;;
	esac
    fi
    o=${o}${sub}/${tmp}
    mkdir -p `dirname $o`
    echo $o
}

# --- Download a single file -----------------------------------------
#  $1: Source file
#  $2: Output directory
#  $3: Number
#  $4: Current count
#  $5: Maximum
#  $6: redirect
#  $7: noact flag
#  $8: multi file flag 
_download_file()
{
    local source=$1 ; shift 
    local store=$1  ; shift 
    local r=$1      ; shift
    local cur=$1    ; shift
    local max=$1    ; shift
    local redir=$1  ; shift
    local noact=$1  ; shift
    local multi=$1  ; shift 

    local o=`_format_out_name "$source" "$store" "$r" "$cur" "$max" "$multi"`
    printf "%4d/%4d: %20s -> %20s ..." $cur $max $source $o 

    mess 2 -n "$source -> $o ... "
    if test -f $o ; then 
	printf "exists\n"
	# mess 2 "exists" 
	# sleep 1
    else
	printf "copying\n"
	# mess 2 -n "copying ... " 
	if test x$noact == x || test $noact -lt 1 ; then 
	    alien_cp alien:${source} file:${o} >> ${redir} 2>&1 
	    fix_perm $o 
	else 
	    sleep 1
	fi
	mess 2 "done"
    fi
    if test x$noact == x || test $noact -gt 0 ; then return 0 ; fi 
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

check_file()
{
    return 0
}
# --- Submit run analysis to background ------------------------------
# $1: Output directory
# $2: Sarting off-set
# $3: Maximum jobs
# $4: Maximum number of files
# $5: noact 
submit_jobs()
{
    mess 3 "submit_jobs out=$1 sta=$2 max=$3 maxf=$4 noact=$5"
    local out=$1 ; shift
    local sta=$1 ; shift 
    local max=$1 ; shift
    local maxf=$1 ; shift
    local noact=$1 ; shift
    local prx=`echo "$spath" | sed 's,//,/,g'` 
    local joblist=
    local counter=0
    mess 5 "Submitting $maxjobs jobs from $sta/$maxf" 
    for i in $@ ; do 
	let cur=$sta+$counter

	local b=`echo $i | sed -e "s,${prx}/*,,"`
	local r=
	mess 10 "submit i=$i b=$b prx=${prx}" 
	case $b in
	    ESDs/*) r=`echo $b | sed -e 's,.*/\([0-9]*\)...\..*/.*,\1,'` ;;
	    *)
		sub=`dirname $b`
		mess 10 "b=$b sub=$sub"
		if test "x$sub" != "x" ; then
		    case $sub in
			[0-9]*/*) sub=`dirname $sub` ;; 
			*/*) sub=`basename $sub`;;
			*) ;;
		    esac
		    mess 10 "b=$b sub=$sub"
		    r=`basename $sub | sed -e "s,/.*,," | sed 's/0*//'`
		fi 
		;;
	esac
	
	let counter=$counter+1

	mess 3 "i=$i out=$out r=$r cur=$cur max=$max noact=$noact " \
	     "(b=$b spath=$spath)"
	download_file "$i" "$out" "$r" "$cur" "$maxf" "$noact" &
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
