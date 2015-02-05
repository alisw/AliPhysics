#!/bin/bash
# 
# BEGIN_MANUAL
# 	Script to help do PWGLF-Forward analsysis using ProofLite
#       =========================================================
# First, one need to figure out what to analyse.  We assume we have
# the ESDs from a real run in some directory - possibly in
# sub-directories, and similar for the MC data. 
# 
# Then, one needs to run this script in set-up mode e.g., 
# 
#   $0 --what=setup  \
#       --name=LHC10c \
#       --run=118560 \
#       --real-dir=/data/alice/data/pp/lhc10c/000118560/pass3 \
#       --real-pattern=AliESDs_*.root \
#       --mc-dir=/data/alice/data/pp/lhc10c/sim/lhc13d4/118560 \
#       --mc-pattern=root_archive.zip@AliESDs.root 
# 
# Note, all the settings are written to the file .config in the
# current directory, so you do not need to give the parameters at
# subsequent steps.  Note, you need a valid AliEn token to at this
# point to get the acceptance corrections.  The run number specified
# is only used for getting the acceptance correction.
# 
# Note, the use of the ZIP archives root_archive.zip and the sub-part
# specification @AliESDs.root for MC data. 
# 
# Next, we need to generate the corrections.  Do 
# 
#   $0 --what=corr 
# 
# and wait for the jobs to finish and terminate. Next, we need to
# extract and upload the corrections to our local corrections folder
# 
#   $0 --what=corr --step=upload 
#   $0 --what=corr --step=draw
# 
# If you already have the corrections, you can pass the option
# --corrections in the setup phase and skip this step.
# 
# Now we can submit our AOD generation jobs.  Do 
# 
#   $0 --what=aod 
# 
# and wait for the jobs to finish and terminate.  If you need to pass
# additional options to the train, one can do so after the special
# option -- e.g., to limit the number of events to 100000, do
# 
#   $0 --what=aod -- --events=100000
# 
# Next, we need to draw the summary results
# 
#   $0 --what=aod --step=draw 
# 
# Now, we should do the dN/deta analysis.  Do 
# 
#   $0 --what=dndeta
# 
# and wait for the jobs to finish and terminate.  Again, additional
# options to the train can be passed after --.  If you passed the
# option --sys=1 in the setup phase, then this will run 3 jobs for
# real and MC each - one for INEL, INEL>0, and NSD (V0-AND).  Next, we
# need to draw the summary and final plot
# 
#   $0 --what=dndeta --step=draw 
# 
# To generate the P(Nch) data, do 
# 
#   $0 --what=multdists
#   $0 --what=multdists --step=draw 
# 
# To collect all PDFs into a single directory do 
# 
#   $0 --what=collect
# 
# Enjoy.
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

run=
real_dir=
real_pat=
real_idx=
mc_dir=
mc_pat=
mc_idx=
my_real_dir=
my_mc_dir=
par=0

# === Various functions ==============================================


# === Implement base functions =======================================
# --- Usage ----------------------------------------------------------
setup_usage()
{
    cat <<EOF
  -r,--run=NUMBER           Specify run number ($run)
  -d,--real-dir=DIR         Directory holding real data ($real_dir)
  -p,--real-pattern=PATTERN Glob pattern to match when searching ($real_pat)
  -D,--mc-dir=DIR           Directory holding MC data ($mc_dir)
  -P,--mc-pattern=PATTERN   Glob pattern to match when searching ($mc_pat)
  -W,--workers=N            Number of workers ($nwrks)
  -a,--par                  Use par files ($par)
EOF
}
# --- handle setup options -------------------------------------------
handle_setup_option()
{
    local arg="$1" 
    local opt="$2"
    case $arg in 
	-r|--run)          run=$opt      ;; 
	-d|--real-dir)     real_dir=$opt ;;
	-p|--real-pattern) real_pat=$opt ;; 
	-D|--mc-dir)       mc_dir=$opt   ;; 
	-P|--mc-pattern)   mc_pat=$opt   ;; 
	-W|--workers)      nwrks=${opt}  ;;
	-a|--par)          par=1         ;;
	*) echo "$0: [SETUP] Unknown option $arg"  ; exit 1 ;; 
    esac
}

# --- Generic draw ---------------------------------------------------
draw()
{
    _draw $@ 
}
# --- Draw dN/deta ---------------------------------------------------
dndeta_draw()
{
    _dndeta_draw $@ 
}
# --- Extract and upload ---------------------------------------------
extract_upload()
{
    echo "=== Download, extract, and uploade in `basename $PWD` ==="
    extract
    upload
}

# --- Get the grid home dir ------------------------------------------
outputs()
{
    l=`pwd`
    my_real_dir="$l/${name}_aod_${now}"
    my_mc_dir="$l/${name}_mcaod_${now}"
    real_idx="$l/${name}_index_${now}.root"
    mc_idx="$l/${name}_mcindex_${now}.root"
}
# --- Dump the setup -------------------------------------------------
dump_setup()
{
    local out=$1 
    cat >> ${out} <<-EOF 
	# Run analysed 
	run=${run}	
	# Real data 
	real_dir=${real_dir}
	real_pat=${real_pat}
	real_idx=${real_idx}
	mc_dir=${mc_dir}
	mc_pat=${mc_pat}
	mc_idx=${mc_idx}
	# Output 
	my_real_dir=${my_real_dir}
	my_mc_dir=${my_mc_dir}
	par=${par}
	EOF
}

run_for_acc()
{
    if test x$run = "x" || test $run -lt 1; then 
	echo "No run for acceptance correction specified" > /dev/stderr 
	exit 1
    fi
    echo $run
}
# --- Run set-ups ----------------------------------------------------
setup()
{
    echo "Calling _setup with $@"
    _setup $@ 

    # create index - unless there's one in the input directory - then 
    # take that and copy here 
    # [We'd like to link only, but ChainBuilder needs to be updated for that] 
    if test "x$real_dir" != x && test -f ${real_dir}/index.root ; then 
	rm -f ${real_idx}
	cp ${real_dir}/index.root ${real_idx} 
    else
	if test ! -f ${real_idx} ; then 
	    index ${real_idx} ${real_dir} "${real_pat}" 0
	fi
    fi
    
    # create index - unless there's one in the input directory - then 
    # take that and copy here 
    # [We'd like to link only, but ChainBuilder needs to be updated for that] 
    if test "x$real_dir" != x && test -f ${mc_dir}/index.root ; then 
	rm -f ${mc_idx}
	cp ${mc_dir}/index.root ${mc_idx}
    else
	if test ! -f ${mc_idx} ; then 
	    index ${mc_idx} ${mc_dir} "${mc_pat}" 0
	fi
    fi
}    

# --- Check settings -------------------------------------------------
check_setup()
{
    if test "x$run" = "x" || test $run -lt 1 ; then 
	echo "Run not specified, or invalid ($run)" > /dev/stderr 
	exit 1
    fi
    if test "x$real_dir" = "x" ; then 
	echo "No real data directory specified" > /dev/stderr 
	exit 1
    fi
    if test "x$mc_dir" = "x" ; then 
	echo "No MC data directory specified" > /dev/stderr 
	# exit 1
    fi
    if test "x$real_pat" = "x" ; then 
	echo "No real data pattern specified" > /dev/stderr 
	exit 1
    fi
    if test "x$mc_pat" = "x" ; then 
	echo "No MC data pattern specified" > /dev/stderr 
	# exit 1
    fi
    
    ncpu=`cat /proc/cpuinfo|sed -n 's/^processor[ \t]*: \(.*\)/\1/p'|wc -l`
    if test "x$nwrks" = "x" || test $nwrks -lt 1 ; then 
	let nwrks=7*$ncpu/10
	echo "Setting number of workers to $nwrks / $ncpu"
    fi
}
# --- Run clean-up ----------------------------------------------------
cleanup()
{
    rm -f ${real_idx} ${mc_idx} 
    _cleanup
}


# === Script specific functions ======================================
# --- Create the index -----------------------------------------------
index()
{
    local o=$1 ; shift
    local d=$1 ; shift 
    local p=$1 ; shift 
    local m=$1 ; shift 

    t=real
    if test $m -gt 0 ; then 
	t=MC
    fi
    if test "x$d" = "x" ;then 
	echo "No input specified for index for $t data" > /dev/stderr 
	return;
    fi
    if test ! -d $d && test ! -f $d ;then 
	echo "Specified input ($d) for $t data is not a directory or file " \
	    > /dev/stderr 
	return;
    fi
    if test $m -gt 0 ; then 
	n="&mc"
    fi
    aliroot -l -b <<EOF
.L $ALICE_PHYSICS/PWGLF/FORWARD/trains/ChainBuilder.C++
ChainBuilder::CreateCollection("${o}", "file://${d}?recursive&scan&pattern=${p}${n}#esdTree");
.q
EOF
    
}
 

# --- Show the setup -------------------------------------------------
print_setup()
{
    cat <<-EOF
	Run:			${run}
	Real data:
	  Directory:	        ${real_dir}
	  Pattern:	        ${real_pat}
	  Output:		${my_real_dir}
	MC data:
	  Directory:		${mc_dir}
	  Pattern:		${mc_pat}
	  Output:		${my_mc_dir}
	Use PAR files:		${par}
	Number of workers:      ${nwrks}/${ncpu}
	EOF
}

# --- Make URI -------------------------------------------------------
url_opts()
{
    local mc=$1 	; shift
    local type=$1 	; shift 
    local trig=$1 	; shift
    local uopt="mode=default&workers=${nwrks}"
    local tree=esdTree 

    local inp=${real_idx}
    if test $mc -gt 0 ; then 
	inp=${mc_idx}
    fi

    case $type in 
	*dndeta|*multdists)
	    tree=aodTree 
	    # Modify for input dir for our files
	    inp=$my_real_dir/AliAOD.root
	    if test $mc -gt 0 ; then 
		inp=$my_mc_dir/AliAOD.root
	    fi
	    ;;
    esac
    if test ! -f $inp ; then 
	echo "No input for $nme, giving up" > /dev/stderr 
	return 
    fi
    if test $mc -gt 0; then 
	uopt="${uopt}&mc"
    fi
    if test $par -gt 0 ; then 
	uopt="${uopt}&par=tasks"
    fi
    if test x$uuopts != x ; then 
	uopt="${uopt}&${uuopts}"
    fi
    # PROOF-lite URL form:
    # 
    #  lite://<datadir_or_list>[?<options>][#<treeName]
    # 
    # Options:
    #  clear=PKGS                 Clear packages ','-separated
    #  mc                         Assume simulation input
    #  mode=default|rec|sim       AliROOT mode 
    #  par=tasks|all              Use par files 
    #  pattern=GLOB               File name pattern 
    #  recursive                  Recursive scan [true]
    #  reset=soft|hard            Reset cluster [hard]
    #  workers=N[x]               Number of workers to use [8]
    #  wrapper=CMD                Wrapper command []

    url="lite://${inp}?${uopt}#${tree}"

}


# === Procedual code =================================================
runIt $@

#
# EOF
#

