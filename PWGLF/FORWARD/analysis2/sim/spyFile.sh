#!/bin/bash

pid=
file=
out=
tee=0
twait=30


# --- Call on errors -------------------------------------------------
error()
{
    echo "$@" > /dev/stderr 
    exit 1
}

# --- Help -----------------------------------------------------------
usage() 
{
    cat <<EOF
Usage: $0 -p PID -f FILE [OPTIONS]

Options:
	-p,--pid	PID		Process identifier 
	-f,--file	FILENAME	File to monitor 
	-o,--out	NAME		Name of stored file 
	-w,--wait	SECONDDS	How many seconds to wait
	-P,--pipe			Use a tee when writing to file 
	-h,--help			This help

Output file is PID_NAME.log.  If NAME is not specified it defaults to
the base name of FILENAME minus possible ending .log. 

EOF
} 

# --- Process command line ------------------------------------------- 
while test $# -gt 0 ; do 
    case $1 in 
	-p|--pid)	pid=$2		; shift	;;
	-f|--file)	file=$2		; shift ;;
	-o|--out)	out=$2		; shift ;; 
	-w|--wait)	twait=$2	; shift ;;
	-P|--pipe)	tee=1		;;
	-h|--help)	usage		; exit 0 ;; 
	*)	error "Unknown option: $1"	;;
    esac
    shift
done

# --- Check settings -------------------------------------------------
if test "X$pid"  = "X" ; then error "No PID specified" ; fi
if test "X$file" = "X" ; then error "No file specified" ; fi
if test "X$out"  = "X" ; then 
    out=`basename $file .log` 
fi
	    
# --- Make derived variables -----------------------------------------
tmp=${pid}_${out}.tmp 
log=${pid}_${out}.log

# --- Check grid authenticity ----------------------------------------
if test ! -f /tmp/gclient_env_${UID} ; then 
    error "No AliEn environment file"
fi
. /tmp/gclient_env_${UID}
alien-token-info | grep -q "Token is still valid"
if test $? -ne 0 ; then 
    echo "Token not valid, please re-new" > /dev/stderr 
    exit 1
fi


# --- Main loop (eternal) --------------------------------------------
while true ; do 
    if test $tee -gt 0 ; then 
	echo "=== Executing alien_spy $pid $file | tee $tmp ..."
	alien_spy $pid $file | tee $tmp
    else
	echo "=== Executing alien_spy $pid $file > $tmp ..."
	alien_spy $pid $file > $tmp
    fi
    l1=`stat -c %s $tmp`
    l2=0
    if test -f $log ; then 
	l2=`stat -c %s $log` 
    fi
    if test $l1 -ge $l2 ; then 
	mv $tmp $log
    else 
	echo "New download smaller than old - exiting"
	break 
    fi 
    sleep $twait
done


#
# EOF
#
