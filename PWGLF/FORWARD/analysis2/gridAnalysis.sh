#!/bin/bash
# 
# BEGIN_MANUAL
# 	Script to help do PWGLF-Forward analsysis using AliEn
#       =====================================================
# First, one need to figure out what to analyse.  Visit the MonAlisa
# pages and find
# 
#   * The list of runs and put that in a file - say runs.list
#   * The directory where the ESD files are stored 
#     - and the pattern that will match these for all runs 
#   * The directory where the MC ESD files are stored 
#     - and the pattern that will match these for all runs 
# 
# Then, one needs to run this script in set-up mode e.g., 
# 
#   $0 --what=setup  \
#       --name=LHC10h \
#       --real-dir=/alice/data/2010/LHC10h \
#       --real-pattern=ESDs/pass2/*/AliESDs.root \
#       --mc-dir=/alice/sim/LHC10h11 \
#       --mc-pattern=*/AliESDs.root \
#       --runs=LHC10.list \
#       --par
# 
# Note, all the settings are written to the file .config in the
# current directory, so you do not need to give the parameters at
# subsequent steps.  As an alternative to giving the parameters, one
# can create the file by hand.
# 
# Next, we need to generate the corrections.  Do 
# 
#   $0 --what=corrs 
# 
# and wait for the jobs to finish and terminate.  If 'watching' is
# turned off, one can also monitor output on MonAlisa, and then when
# enough has finished, execute
# 
#   $0 --what=corrs --step=terminate 
# 
# enough times to get the final merged result.  Next, we need to
# extract and upload the corrections to our local corrections folder
# 
#   $0 --what=corrs --step=upload 
#   $0 --what=corrs --step=draw
# 
# Now we can submit our AOD generation jobs.  Do 
# 
#   $0 --what=aod
# 
# and wait for the jobs to finish and terminate.  If 'watching' is
# turned off, one can also monitor output on MonAlisa, and then when
# enough has finished, execute
# 
#   $0 --what=aod --step=terminate 
# 
# enough times to get the final merged result.  Next, we need to
# download the results and we draw the summary results
# 
#   $0 --what aod --step=draw 
# 
# Now, we should do the dN/deta analysis.  Do 
# 
#   $0 --what=dndeta
# 
# and wait for the jobs to finish and terminate.  If 'watching' is
# turned off, one can also monitor output on MonAlisa, and then when
# enough has finished, execute
# 
#   $0 --what=dndeta --step=terminate 
# 
# enough times to get the final merged result.  If you passed the
# option --sys=1 in the setup phase, then this will run 3 jobs for
# real and MC each - one for INEL, INEL>0, and NSD (V0-AND). Next, we
# need to download the results and we can draw the summary and final
# plot
# 
#   $0 --what=dndeta --step=draw 
# 
# When running the trains, one can pass additional options to the
# train after the special option -- e.g., 
#   
#   $0 --what=aod -- --verbose=2 --branches --satellite
# 
# Enjoy
# 
# Comments, questions, bugs, flames, suggestions, etc. should be sent
# to Christian Holm Christensen <cholm@nbi.dk>
# 
# END_MANUAL
if test ! -f $ALICE_PHYSICS/PWGLF/FORWARD/analysis2/baseAnalysis.sh  ; then 
    echo "baseAnalysis not found!" > /dev/stderr 
    exit 1
fi
. $ALICE_PHYSICS/PWGLF/FORWARD/analysis2/baseAnalysis.sh 

runs=
mcruns=
par=0
real_dir=
real_pat=
mc_dir=
mc_pat=
my_real_dir=
my_mc_dir=
watch=0

# === Various functions ==============================================


# === Implement base functions =======================================
# --- Usage ----------------------------------------------------------
setup_usage()
{
    cat <<EOF
  -r,--runs=FILENAME        Specify list of runs file ($runs)
  -R,--mc-runs=FILENAME     Specify list of MC runs file ($mcruns)
  -d,--real-dir=ALIEN_DIR   Directory holding real data ($real_dir)
  -p,--real-pattern=PATTERN Glob pattern to match when searching ($real_pat)
  -D,--mc-dir=ALIEN_DIR     Directory holding MC data ($mc_dir)
  -P,--mc-pattern=PATTERN   Glob pattern to match when searching ($mc_pat)
  -W,--watch                Watch for job status and terminate automatically
  -a,--par                  Use par files ($par)
EOF
}
handle_setup_option()
{
    local arg="$1"
    local opt="$2"
    # echo "Grid: Processing '$arg' ('$opt')"
    case $arg in 
	-r|--runs)         runs=$opt     ;; 
	-R|--mc-runs)      mcruns=$opt   ;;
	-d|--real-dir)     real_dir=$opt ;;
	-p|--real-pattern) real_pat=$opt ;; 
	-D|--mc-dir)       mc_dir=$opt   ;; 
	-P|--mc-pattern)   mc_pat=$opt   ;; 
	-W|--watch)        let watch=\!$watch ;;
	-a|--par)          par=1         ;;
	*) echo "$0: [SETUP] Unknown option $arg"  ; exit 1 ;; 
    esac
}    
# --- Extract corrections --------------------------------------------
terminate()
{
    script Terminate.C 
}

# --- Extract and upload ---------------------------------------------
extract_upload()
{
    echo "=== Download, extract, and uploade in `basename $PWD` ==="
    download
    for i in *.zip ; do 
	if test ! -f $i ; then continue ; fi 
	if test $noact -gt 0 ; then continue ; fi
	local d=`basename $i .zip` 
	if test ! -d $d ; then 
	    mkdir $d 
	fi
	cd $d 
	if test ! -f .zip ; then 
	    echo "= Unpacking ../$i"
	    unzip ../$i > /dev/null 2>&1
	    touch .zip 
	fi
	_extract ../Extract.C 
	upload
	cd ..
    done
}

# --- Draw -----------------------------------------------------------
draw()
{
    local dd=`pwd` 
    dd=`basename $dd` 
    echo "=== Drawing in $dd using $1"
    local scr=$1  ; shift    
    case x$scr in 
	x/*) ;; 
	x*)  scr=../$scr ;; 
    esac
    # local args=$1 ; shift
    download 
    for i in *.zip ; do 
	if test "X$i" = "X*.zip" ; then continue ; fi
	echo "--- Will extract $i in $dd"
	d=`basename $i .zip` 
	if test ! -d $d ; then 
	    mkdir -p $d 
	    unzip $i -d $d
	fi
	(cd $d && _draw $scr $@)
    done
}
dndeta_draw()
{
    echo "=== Drawing dN/deta in $1"
    local top=$1 ; shift 
    cd $top
    download 
    for i in *.zip ; do 
	if test "X$i" = "X*.zip" ; then continue ; fi
	echo "--- Will extract $i"
	d=`basename $i .zip` 
	if test ! -d $d ; then 
	    mkdir -p $d 
	    unzip $i -d $d
	fi
	(cd $d && \
	    script ${fwd_dir}/DrawdNdetaSummary.C && \
	    script ../Draw.C)
    done
    cd ..
}

# === Script specific functions ======================================
# --- Extract corrections --------------------------------------------
download()
{
    echo "=== Executing download in `pwd`"
    if test -f .download ; then 
	echo "--- Already downloaded in `basename $PWD`"
	return 0 
    fi
    script Download.C 
    touch .download
}
    
# --- Get the grid home dir ------------------------------------------
outputs()
{
    l=`aliroot -l -b <<EOF
gSystem->RedirectOutput("/dev/null");
TGrid::Connect("alien://");
gSystem->RedirectOutput(0);
std::cout << gGrid->GetHomeDirectory() << std::endl;
EOF`
    my_real_dir="$l/${name}_aod_${now}/output"
    my_mc_dir="$l/${name}_mcaod_${now}/output"
}


run_for_acc()
{
    local r=`grep -v ^# ../$runs | awk '{FS=" \n\t"}{printf "%d\n", $1}' | head -n 1` 
    if test x$r = "x" || test $r -lt 1; then 
	echo "No run for acceptance correction specified" > /dev/stderr 
	exit 1
    fi
    echo $r 
}

# --- Dump the setup -------------------------------------------------
dump_setup()
{
    local out=$1 
    cat >> ${out} <<-EOF 
	# Real data 
	runs=${runs}
	real_dir=${real_dir}
	real_pat=${real_pat}
	# Simulated data
	mcruns=${mcruns}
	mc_dir=${mc_dir}
	mc_pat=${mc_pat}
	# Output directories 
	my_real_dir=${my_real_dir}
	my_mc_dir=${my_mc_dir}
	# Other options 
	par=${par}
	watch=${watch}
	EOF
}

# --- Run set-ups ----------------------------------------------------
setup()
{
    _setup $@
}    


# --- Check settings -------------------------------------------------
check_setup()
{
    check_token

    if test "x$runs" = "x" || test ! -f $runs ; then 
	echo "List of run file $runs not found" > /dev/stderr 
	exit 1
    fi
    if test "x$mcruns" = "x" ; then mcruns=$runs ; fi 
    if test ! -f $mcruns ; then 
	echo "List of MC runs file $mcruns not found" > /dev/stderr 
	exit 1
    fi
    if test "x$real_dir" = "x" ; then 
	echo "No real data directory specified" > /dev/stderr 
	exit 1
    fi
    if test "x$mc_dir" = "x" ; then 
	echo "No MC data directory specified" > /dev/stderr 
	exit 1
    fi
    if test "x$real_pat" = "x" ; then 
	echo "No real data pattern specified" > /dev/stderr 
	exit 1
    fi
    if test "x$mc_pat" = "x" ; then 
	echo "No MC data pattern specified" > /dev/stderr 
	exit 1
    fi

}

# --- Show the setup -------------------------------------------------
print_setup()
{
    cat <<-EOF
	Real data:
	  Run file:	        ${runs}
	  Directory:		${real_dir}
	  Pattern:		${real_pat}
	  Output:		${my_real_dir}
	MC data:
	  Run file:       	${mcruns}
	  Directory:		${mc_dir}
	  Pattern:		${mc_pat}
	  Output:		${my_mc_dir}
	Use PAR files:		${par}
	EOF
}



# --- Make URI -------------------------------------------------------
# Must modify URL 
url_opts()
{
    local mc=$1 	; shift
    local type=$1 	; shift 
    local trig=$1 	; shift

    local uopt="&merge=50&split=50&aliroot=last,regular"

    local dir=$real_dir
    local pat=$real_pat
    local rl=$runs
    local tree=esdTree

    if test $mc -gt 0 ; then 
	dir=$mc_dir
	pat=$mc_pat
	rl=$mcruns
    fi
    case $type in 
	*dndeta|*multdists) 
	    uopt="${uopt}&concat"
	    tree=aodTree 
	    # Modify for input dir for our files
	    dir=$my_real_dir
	    pat="*/AliAOD.root"
	    if test $mc -gt 0 ; then 
		dir=$my_mc_dir
	    fi
	    ;;
    esac
    if test $mc -gt 0; then 
	uopt="${uopt}&mc"
    fi
    if test $par -gt 0 ; then 
	uopt="${uopt}&par"
    fi
    if test x$uuopts != x ; then 
	uopt="${uopt}&${uuopts}"
    fi
    url="alien://${dir}?run=${rl}&pattern=${pat}${uopt}#${tree}"    
}

    
# --- Run the train --------------------------------------------------
# Usage:
# 
allAboard()
{
    local type=$1 ; shift
    local trig=$1 ; shift
    _allAboard "$type" "$trig" --batch $@ 

    if test $watch -lt 1 ; then 
	cat <<-EOF
	Check https://alimonitor.cern.ch/users/jobs.jsp for progress 
	
	Remember to do 
	
	  $0 --what=... --step=terminate
	
	until the final merge stage to get the results. 
	EOF
	
	case $type in 
	    *corr|*esd) 
		cat <<-EOF
	Then, do 
	
	  $0 --what=... --step=upload 
	
	to upload the results to our local corrections store. 
	EOF
		;; 
	    *aod|*dndeta)
		cat <<-EOF
	Then, do 
	 
	  $0 --what=... --step=draw
	
	to get a PDF of the diagnostics histograms and the final plots. 
	EOF
		;;
	esac
    else 
	echo "Now waiting for jobs to finish"
	(cd ${nme}_${now} && \
	    nice aliroot -l -b -x -q Watch.C\(1\) 2>&1 | \
	    tee watch.log > /dev/null &)
    fi
}

# --- Collect a directory --------------------------------------------
collect_files()
{
    local dir=$1 ; shift
    local d=$1 ; shift 
    local M=$1 ; shift 
    local out=$1 ; shift
    local r=$1 ; shift
    local files="$1"

    for ad in $dir/root_archive_* ; do 
	if test ! -d $ad ; then continue ; fi 
	local r=`basename $ad | sed 's/root_archive_0*//'` 
	_collect_files "$ad" "$d" "$M" "$out" "$r" "$files"
    done # for ad in ...
}

    
# === Procedual code =================================================
runIt $@

#
# EOF
#

