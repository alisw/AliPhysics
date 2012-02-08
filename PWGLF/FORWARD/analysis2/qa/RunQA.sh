#!/bin/bash
#
# This is deprecated.  Use RunQAMT.sh instead
#
jobid=0
top=.
verb=0
nodw=0
notr=0
norn=0
maxf=-1
last="unknown"
lock=


export PATH=$PATH:$ALICE_ROOT/PWGLF/FORWARD/analysis2/qa:../scripts

# --- Handle exit ----------------------------------------------------
handle_exit()
{
    echo "Removing $lock"
    rm -f $lock
}
trap handle_exit EXIT 

# --- Usage information ----------------------------------------------
usage()
{ 
cat <<EOF
Usage: $0 -j JOBID [OPTIONS]

Options:
	-h.--help		   This help
	-j,--jobid	JOBID	   Job id from MonALisa                [$jobid]
	-d,--nodownload            Do not download                     [$nodw]
	-t,--notrend               Do not make new trend tree          [$notr]
	-r,--no-run                Do not make per-run info            [$norn]
	-o,--output	DIRECTORY  Where to store the result           [$top]
	-v,--verbose		   Be verbose                          [$verb]
	-m,--max	NUMBER	   Maximum number of runs to get       [$maxf]
	-s,--skip-lines NUMBER	   Number of lines to skip in job list [$skip]
EOF
}     

# --- Messages -------------------------------------------------------
mess()
{
    if test $verb -lt 1 ; then return ; fi 
    last="$@"
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
    Y=`echo $p | sed 's/LHC\(..\).*/\1/'` 
    L=`echo $p | sed "s/LHC${Y}\(.\).*/\1/"` 
    S=`echo $p | sed "s/LHC${Y}${L}//"` 

    dprod=LHC${Y}${L}
    case x$S in 
	x) ;; 
	x_*) ;; 
	*) dprod=${dprod}${S} ;; 
    esac
    
    case x$r in 
	xQA*) q=`echo $r | sed 's/QA//'` ;; 
	x) ;;
	*) ;;
    esac
	    
    opts="-p $p -P $P"
    if test "x$R" != "x" ; then opts="$opts -R $R" ; fi 
    if test "x$q" != "x" ; then opts="$opts -q $q" ; fi 
    if test "x$Q" != "x" ; then opts="$opts -Q $Q" ; fi 

    dest="${dprod}/${Q}pass${P}"
    if test "x$R" != "x" ; then dest="${dest}_${R}" ; fi 

}

# --- Get the options ------------------------------------------------
skip=1
jobUrl="http://alimonitor.cern.ch/prod/jobs.jsp?t="
get_opts() { 
    wget -q ${jobUrl}${jobid} -O job.html
    p=`grep "/catalogue/index.jsp?path" job.html | head -n $skip | tail -n 1 | sed -e 's,.*/alice/data/\([^<]*\)<.*,\1,' | tr '/' ' '` 
    rm -f job.html 

    get_parts $p
} 

# --- Trapping errors ------------------------------------------------
handle_err()
{
    echo "Got an error: $last"
    exit 1
} 
enable_trap ()
{
    # echo "Enabling trapping errors"
    trap handle_err ERR 
    # trap -p ERR
}
disable_trap ()
{
    # echo "Disabling trapping errors"
    trap - ERR
    # trap -p ERR
}


# --- Deprecated -----------------------------------------------------
cat <<EOF
This script is deprecated.  You should use the more advanced script RunQAMT.sh 
which process the runs in parallel and in general is more flexible
EOF
read -n 1 -i n -p "Do you want to continue [yN]? "
case $REPLY in 
    y|Y) echo "" ;; 
    *) exit 0 ;; 
esac


# --- comamnd line ---------------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)	  usage   ; exit 0 ;;
	-j|--jobid)	  jobid=$2; shift  ;;
	-o|--output)	  top=$2  ; shift  ;; 
	-d|--no-download) nodw=1  ;;
	-t|--no-trend)    notr=1  ;;
	-r|--no-run)      norn=1  ;; 
	-v|--verbose)	  verb=1  ;;
	-s|--skip-lines)  skip=$2 ; shift ;; 
	-m|--max)         maxf=$2 ; shift ;; 
	*) echo "Unknown option $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift 
done

# --- Sanity ---------------------------------------------------------
if test $jobid -lt 1 ; then 
    echo "No JOBID specified" > /dev/stderr 
    exit 1
fi

# --- Extract options -----------------------------------------------
get_opts $jobid
if test "x$opts" = "x" ; then 
    echo "No options found" 
    rm -f $lock
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
Skip lines:	  $skip
Max # of runs:    $maxf
----------------------------------------------------------------------
EOF

# --- Lock -----------------------------------------------------------
lock=${top}/${dest}/.lock
if test -f $lock ; then 
    echo "Another QA process is already running:" > /dev/stderr 
    echo "Content of ${top}/${dest}/.lock:" > /dev/stderr 
    cat $lock > /dev/stderr 
    trap - EXIT
    exit 1
fi

now=`date`
cat <<EOF > $lock
Process: $$ 
User:    $USER 
Start:   $now 
EOF


# --- Download the files ---------------------------------------------
if test $nodw -lt 1 ; then 
    mess "Running getQAResults.sh $opts -d $top -n "
    getQAResults.sh $opts -d $top -T -v -v -i -m $maxf
else 
    mess "Not downloading"
fi

# --- Now run the QA code -------------------------------------------

mess "Now running code"
enable_trap
savdir=`pwd`
mess "Change directory to $top/$dest"
cd $top/${dest} 

disable_trap
trap - ERR
trap -p ERR 
mess "List of trend_<x>.html files"
idx=`ls trend_*_*.html 2>/dev/null`
mess "Removing indeces"  
rm -f index.html
for i in $idx ; do 
    if test -f $i ; then rm -f $i ; fi
done 

what=3
if test $notr -gt 0 ; then let what=$what^0x2 ; fi
if test $norn -gt 0 ; then let what=$what^0x1 ; fi

scr=$ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/RunQA.C
mess "Running root -l -b -q ${scr}\(\".\",1,-1,$what\)"

enable_trap
root -l -b -q ${scr}\(\".\",1,-1,$what\)
trap - ERR
disable_trap

idx=`ls trend_*_*.html 2>/dev/null` 
rm -f index.html
for i in $idx ; do 
    echo "making index.html point to $i"
    cat $i | \
	sed -e "s,index.html,../index.html," \
	    -e "s,!--JOBID--,a target='_blank' href='${jobUrl}${jobid}'>Job</a," \
	> index.html  
    cp index.html $i 
done 
if test ! -f index.html ; then 
    echo "No index file found" 
    exit 1
fi
chmod g+rw  index.html
chmod g+rw .  > /dev/null 2>&1

style= 
if test -f $ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/style.css ; then 
    style=$ALICE_ROOT/PWGLF/FORWARD/analysis2/qa/style.css
elif test -f $savdir/style.css ; then 
    style=$savdir/style.css 
fi 
   
if test x$style != x ; then 
    rm -f style.css 
    cp $style .
    chmod g+rw style.css
fi

cat <<EOF

Finished running QA for jobid $jobid.  
Output is stored in  $top/${dest} 

EOF

# --- Make index.html ------------------------------------------------
cd ..
date=`date`
period=`pwd`
period=`basename $period`
rm -f index.html
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
<div class='jobid'>
   
</div>
<div class='back'>
  <a href="../">Back</a>
</div>
<div class='change'>
  Last update: $date
</div>
</body>
</html>
EOF
chmod g+rw index.html 
chmod g+rw .  > /dev/null 2>&1
if test x$style != x ; then 
    rm -f style.css
    cp $style .
    chmod g+rw style.css
fi

# --- Make index.html ------------------------------------------------
cd ..
rm -f index.html
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
<div class='change'>
Last update: $date
</div>
</body>
</html>
EOF
chmod g+rw index.html 
chmod g+rw .  > /dev/null 2>&1
if test x$style != x ; then 
    rm -f style.css
    cp $style . 
    chmod g+rw style.css
fi

cd $savdir
rm -f $lock

# 
# EOF
#




