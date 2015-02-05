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
serv=
real_ds=
mc_ds=
my_real_ds=
my_mc_ds=
nwrkrs=0
par=0
batch=0
# === Various functions ==============================================


# === Implement base functions =======================================
# --- Usage ----------------------------------------------------------
setup_usage()
{
    cat <<EOF
  -r,--run=NUMBER           Specify run number ($run)
  -d,--real-ds=NAME         Data set of real data ($real_ds)
  -D,--mc-ds=NAME           Data set of MC data ($mc_dir)
  -P,--server=ADDRESS	    Proof server ($serv)
  -W,--workers=N            Number of workers ($nwrks)
  -a,--par                  Use par files ($par)
  -b,--batch                Do not show GUI
EOF
}
# --- handle setup options -------------------------------------------
handle_setup_option()
{
    local arg="$1" 
    local opt="$2"
    # echo "Handle arg=$arg opt=$opt"
    case $arg in 
	-r|--run)          run=$opt      ;; 
	-d|--real-ds)      real_ds=$opt  ;;
	-D|--mc-ds)        mc_ds=$opt    ;; 
	-W|--workers)      nwrks=${opt}  ;;
	-P|--server)       serv=${opt}   ;;
	-a|--par)          par=1         ;;
	*) echo "$0: [SETUP] Unknown option $arg"  ; exit 1 ;; 
    esac
}

# --- Get the grid home dir ------------------------------------------
outputs()
{
    # We should get the PROOF group and user name here
    # l=`id -n -g`/`id -n -u`
    l=`aliroot -l -b <<EOF
gSystem->RedirectOutput("/dev/null");
TProof::Reset("$serv", false);
TProof::Open("$serv");
gSystem->RedirectOutput(0);
Printf("%s/%s", gProof->GetGroup(), gProof->GetUser());
EOF`
    my_real_ds="/$l/${name}_aod_${now}"
    my_mc_ds="/$l/${name}_mcaod_${now}"
}
# --- Dump the setup -------------------------------------------------
dump_setup()
{
    local out=$1 
    cat >> ${out} <<-EOF 
	# Run analysed 
	run=${run}	
	# Real data 
	real_ds=${real_ds}
	mc_ds=${mc_ds}
	# Output 
	my_real_ds=${my_real_ds}
	my_mc_ds=${my_mc_ds}
	serv=${serv}
	nwrks=${nwrks}
	par=${par}
	batch=${batch}
	EOF
}

# ---- Get run number to use for Acceptance map ----------------------
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
}    

# --- Check settings -------------------------------------------------
check_setup()
{
    if test "x$run" = "x" || test $run -lt 1 ; then 
	echo "Run not specified, or invalid ($run)" > /dev/stderr 
	exit 1
    fi
    if test "x$real_ds" = "x" ; then 
	echo "No real data set specified" > /dev/stderr 
	exit 1
    fi
    if test "x$mc_ds" = "x" ; then 
	echo "No MC data set specified" > /dev/stderr 
	# exit 1
    fi
}



# --- Show the setup -------------------------------------------------
print_setup()
{
    cat <<-EOF
	Run:			${run}
	Real data set:
	  Data set:		${real_ds}
	  Output:		${my_real_ds}
	MC data:
	  Data set:		${mc_ds}
	  Output:		${my_mc_ds}
	Use PAR files:		${par}
	Number of workers:      ${nwrks}
	Server:                 ${serv}
	EOF
}

# --- Make URI -------------------------------------------------------
url_opts()
{
    local mc=$1 	; shift
    local type=$1 	; shift 
    local trig=$1 	; shift
    local uopt="mode=default"
    local tree=esdTree 

    local inp=${real_ds}
    local oup=`basename ${my_real_ds}`
    if test $mc -gt 0 ; then 
	inp=${mc_ds}
	oup=`basename ${my_mc_ds}`
    fi

    case $type in 
	*aod)
	    uopt="${uopt}&dsname=${oup}"
	    ;; 
	*dndeta|*multdists)
	    tree=aodTree 
	    # Modify for input dir for our files
	    inp=$my_real_ds
	    if test $mc -gt 0 ; then 
		inp=$my_mc_ds
	    fi
	    ;;
    esac
    if test x$inp = x; then 
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
    if test "x$nwrks" != "x" && test $nwrks -gt 0 ; then 
	uopt="${uopt}&workers=${nwrks}"
    fi
    case $server in 
	alice-caf.cern.ch) uopts="${uopts}&reset=hard" ;; 
    esac 
	    
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

    url="proof://${serv}/${inp}?${uopt}#${tree}"

}

# --- Run the train --------------------------------------------------
# Usage:
# 
allAboard()
{
    local type=$1 ; shift
    local trig=$1 ; shift
    local lopts=
    if test $batch -gt 0 ; then 
	lopts="--batch"
    fi
    _allAboard "$type" "$trig" $lopts $@ 
}

# === Procedual code =================================================
runIt $@

#
# EOF
#

