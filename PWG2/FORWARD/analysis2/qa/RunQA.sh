#!/bin/bash

jobid=0
top=.
verb=0
nodw=0
notr=0

export PATH=$PATH:$ALICE_ROOT/PWG2/FORWARD/analysis2/qa:../scripts

usage()
{
cat <<EOF
Usage: $0 -j JOBID [OPTIONS]

Options:
	-h.--help		   This help
	-j,--jobid	JOBID	   Job id from MonALisa       [$jobid]
	-d,--nodownload            Do not download            [$nodw]
	-t,--notrend               Do not make new trend tree [$notr]
	-o,--output	DIRECTORY  Where to store the result  [$top]
	-v,--verbose		   Be verbose                 [$verb]
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
    rm -f job.html 
    get_parts $p
}

# --- comamnd line ---------------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)	  usage   ; exit 0 ;;
	-j|--jobid)	  jobid=$2; shift  ;;
	-o|--output)	  top=$2  ; shift  ;; 
	-d|--no-download) nodw=1  ;;
	-t|--no-trend)    notr=1  ;;
	-v|--verbose)	  verb=1  ;;
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

# --- Display job setting --------------------------------------------
cat <<EOF
----------------------------------------------------------------------
Job id:		  $jobid
Top directory:	  $top
Don't download:   $nodw
Don't make tree:  $notr
Download options: $opts
Destination:      $dest
----------------------------------------------------------------------
EOF

# --- Exit on errors -------------------------------------------------
set -e 

# --- Download the files ---------------------------------------------
if test $nodw -lt 1 ; then 
    mess "Running getQAResults.sh $opts -d $top -n "
    getQAResults.sh $opts -d $top -T -v -v 
fi

# --- Now run the QA code -------------------------------------------
savdir=`pwd`
cd $top/${dest} 
what=3
if test $notr -gt 0 ; then what=2 ; fi
scr=$ALICE_ROOT/PWG2/FORWARD/analysis2/qa/RunQA.C
mess "Running root -l -b -q ${scr}\(\".\",1,-1,$what\)"
root -l -b -q ${scr}\(\".\",1,-1,$what\)
idx=`ls trend_*_*.html 2>/dev/null` 
for i in $idx ; do 
    echo "making index.html point to $i"
    ln -fs $i index.html  
done 
chmod g+rw  index.html
if test -f $savdir/style.css ; then 
    cp $savdir/style.css .
fi

cat <<EOF

Finished running QA for jobid $jobid.  
Output is stored in  $top/${dest} 

EOF

# --- Make index.html ------------------------------------------------
cd ..
period=`pwd`
period=`basename $period`
cat <<EOF > index.html
<html>
<head>
<title>$period</title>
<LINK REL="stylesheet" href="style.css">
</head>
<body>
<h1>$period</h1>
<ul>
EOF
for i in * ; do 
    if test ! -d $i ; then continue ; fi 
    echo "<li><a href='$i'>$i</a></li>" >> index.html 
done
cat <<EOF >> index.html
</ul>
<div class='back'><a href="../">Back</a></div>
</body>
</html>
EOF
chmod g+rw index.html 
if test -f $savdir/style.css ; then 
    cp $savdir/style.css .
fi

# --- Make index.html ------------------------------------------------
cd ..
cat <<EOF >index.html
<html>
<head>
<title>QA for the FMD</title>
<LINK REL="stylesheet" href="style.css">
</head>
<body>
<h1>QA for the FMD</h1>
<p>
For more information, see <a href='https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FMDQA'>TWiki pages</a>
</p>
<ul>
EOF
for i in * ; do 
    if test ! -d $i ; then continue ; fi 
    echo "<li><a href='$i'>$i</a></li>" >> index.html 
done
cat <<EOF >> index.html
</ul>
</body>
</html>
EOF
chmod g+rw index.html 
if test -f $savdir/style.css && test `pwd` != $savdir; then 
    cp $savdir/style.css .
fi

cd $savdir

# 
# EOF
#




