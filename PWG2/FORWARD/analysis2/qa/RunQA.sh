#!/bin/bash

jobid=0
top=.
verb=0

export PATH=$PATH:$ALICE_ROOT/PWG2/FORWARD/analysis2/qa:../scripts

usage()
{
cat <<EOF
Usage: $0 -j JOBID [OPTIONS]

Options:
	-h.--help		   This help
	-j,--jobid	JOBID	   Job id from MonALisa
	-o,--output	DIRECTORY  Where to store the result
	-v,--verbose		   Be verbose 
EOF
} 

mess()
{
    if test $verb -lt 1 ; then return ; fi 
    echo $@
}

# --- Get parts of the passed path -----------------------------------
get_parts()
{
    y=$1 ; shift 
    p=$1 ; shift 
    r=$1 ; shift 
    e=$1 ; shift 
    x=$1 ; shift 
    r=$1 
    
    P=`echo $x | sed 's/.*pass\([0-9]*\).*/\1/'`
    R=`echo $x | sed -n "s/.*pass${P}_//p"` 
    Q=`echo $x | sed -n 's/pass.*//p'`
    
    case x$r in 
	xQA*) q=`echo $r | sed 's/QA//'` ;; 
	x) ;;
	*) ;;
    esac
	    
    opts="-p $p -P $P"
    if test "x$R" != "x" ; then opts="$opts -R $R" ; fi 
    if test "x$q" != "x" ; then opts="$opts -q $q" ; fi 
    if test "x$Q" != "x" ; then opts="$opts -Q $Q" ; fi 

    dest="${p}/${Q}pass${P}"
    if test "x$R" != "x" ; then dest="${dest}_${R}" ; fi 

}

# --- Get the options ------------------------------------------------
get_opts() {
    wget -q http://alimonitor.cern.ch/prod/jobs.jsp?t=${jobid} -O job.html
    p=`grep "/catalogue/index.jsp?path" job.html | head -n 1 | sed -e 's,.*/alice/data/\([^<]*\)<.*,\1,' | tr '/' ' '` 
    get_parts $p
}

# --- comamnd line ---------------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)	usage 	; exit 0 ;;
	-j|--jobid)	jobid=$2; shift  ;;
	-o|--output)	top=$2 ; shift  ;; 
	-v|--verbose)	verb=1   ;;
	*) echo "Unknown option $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift 
done

# --- Sanity ---------------------------------------------------------
if test $jobid -lt 1 ; then 
    echo "No JOBID specified" > /dev/stderr 
    exit 1
fi

# --- Extract options ------------------------------------------------
get_opts $jobid
mess "Options for getQAresults: $opts, destingation: $dest" 
if test "x$opts" = "x" ; then 
    echo "No options found" 
    exit 1
fi

set -e 

# --- Download the files ---------------------------------------------
mess "Running getQAResults.sh $opts -d $top -n "
getQAResults.sh $opts -d $top -T -v -v 

# --- Now run the QA code -------------------------------------------
savdir=`pwd`
cd $top/${dest} 
root -l -b -q $ALICE_ROOT/PWG2/FORWARD/analysis2/qa/RunQA.C
cd $savdir

cat <<EOF

Finished running QA for jobid $jobid.  
Output is stored in  $top/${dest} 
EOF

# 
# EOF
#




