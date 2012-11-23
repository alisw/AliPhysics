#!/bin/bash
# 
# BEGIN_MANUAL
# 
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
#       --name=LHC10h --sys=pbpb --snn=2760 --field=5 \
#       --real-dir=/alice/data/2010/LHC10h \
#       --real-pattern=ESDs/pass2/*/AliESDs.root \
#       --mc-dir=/alice/sim/LHC10h11 \
#       --mc-pattern=*/AliESDs.root \
#       --runs=LHC10.list --par
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
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what=corrs --step=terminate 
# 
# enough times to get the final merged result.  Next, we need to
# extract and upload the corrections to our local corrections folder
# 
#   $0 --what=corrs --step=upload 
# 
# Now we can submit our AOD generation jobs.  Do 
# 
#   $0 --what=aods 
# 
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what=aods --step=terminate 
# 
# enough times to get the final merged result.  Next, we need to
# download the results and we draw the summary results
# 
#   $0 --what aods --step=draw 
# 
# Now, we should do the dN/deta analysis.  Do 
# 
#   $0 --what=dndetas
# 
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what=dndetas --step=terminate 
# 
# enough times to get the final merged result.  Next, we need to download 
# the results and we can draw the summary and final plot
# 
#   $0 --what=dndetas --step=draw 
# 
# Enjoy
# 
# END_MANUAL

runs=
name=
sys=1
snn=900
field=5
corrs=corrections
dotconf=.config
here=${PWD}
par=0
noact=0
aliroot="&aliroot=v5-03-75pATF-AN"
# root="root=v5-34-02-1"
fwd_dir=$ALICE_ROOT/PWGLF/FORWARD/analysis2

real_dir=
real_pat=
mc_dir=
mc_pat=
my_real_dir=
my_mc_dir=

# === Various functions ==============================================
# --- Usage ----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 --what OPERATION [OPTIONS]

Options:
  -r,--runs=FILENAME        Specify list of runs file ($runs)
  -n,--name=STRING          Base name of jobs ($name)
  -S,--sys=SYSTEM           Collision system ($sys)
  -E,--snn=ENERGY           Center of mass energy per nuclean pair ($snn)
  -F,--field=FIELD          L3 magnetic field ($field)
  -d,--real-dir=ALIEN_DIR   Directory holding real data ($real_dir)
  -p,--real-pattern=PATTERN Glob pattern to match when searching ($real_pat)
  -D,--mc-dir=ALIEN_DIR     Directory holding MC data ($mc_dir)
  -P,--mc-pattern=PATTERN   Glob pattern to match when searching ($mc_pat)
  -s,--step=STEP            Run stage ($step)
  -w,--what=TRAINS          What to do 
  -a,--par                  Use par files ($par)
  -M,--man                  Show the manual  
  -N,--noact                Show what will be done 

TRAINS is one of

  clean       Clean directory
  setup       Do intial setup 
  corrs       Generate corrections 
  aods        Generate AODs 
  dndeta      Generate dNdeta 

and must be executed in that order.  STEP is one of 

  full        Run the analysis 
  terminate   Terminate the job (may need iterations)
  upload      Upload corrections (only for TRAINS=corrs)
  draw        Draw (partial) results (not for TRAINS=corrs)
EOF
}

# --- Manual ---------------------------------------------------------
manual()
{
    grep ^# $0 | \
	sed -n -e '/BEGIN_MANUAL/,/END_MANUAL/ p' | \
	sed -e '/\(BEGIN\|END\)_MANUAL/ d' -e 's/^# //' -e "s,\$0,$0,"
}

# === Utilities to execute scripts ===================================
# --- Run script -----------------------------------------------------
script()
{
    local scr=$1 ; shift 
    local args=$1 ; shift
    echo "Will run aliroot -l -b -q $scr($args)"
    if test $noact -gt 0 ; then return ; fi
    aliroot -l -b <<EOF
.x $scr($args)
.q
EOF
}
# --- Run acceptance generation --------------------------------------
accGen()
{
    local run=$1 
    script ${fwd_dir}/corrs/ExtractAcceptance.C \
	"${run},${sys},${snn},${field}"
}

# --- Extract corrections --------------------------------------------
terminate()
{
    script Terminate.C 
}
# --- Extract corrections --------------------------------------------
download()
{
    test -f .download && return 0 
    script Download.C 
    touch .download
}
# --- Extract corrections --------------------------------------------
extract()
{
    test -f .extract && return 0 
    script Extract.C 
    touch .extract
}
# --- Upload a file --------------------------------------------------
upload()
{
    test -f .upload && return 0 
    script Upload.C \"file://${here}/${name}_corrs_${now}\"
    touch .upload 
}
# --- Extract and upload ---------------------------------------------
extract_upload()
{
    echo "Download, extract, and uploade in `basename $PWD`"
    download
    for i in *.zip ; do 
	if test ! -f $i ; then continue ; fi 
	echo "Extracting and uploading from $i"
	if test $noact -gt 0 ; then continue ; fi
	rm -rf tmp
	(mkdir -p tmp && \
	    cd tmp && \ 
	    unzip ../$i && \
		script ../Extract.C "" "" && 
	    upload)
    done
}

# --- Draw -----------------------------------------------------------
draw()
{
    local scr=$1 
    download 
    for i in *.zip ; do 
	d=`basename $i .zip` 
	mkdir -p $d 
	unzip $i -d $d
	(cd $d && script $scr)
    done
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
    my_real_dir="$l/${name}_${now}_aod/output"
    my_mc_dir="$l/${name}_${now}_mcaod/output"
}

# === Trains =========================================================
# --- Run set-ups ----------------------------------------------------
setup()
{
    run_for_acc=`cat $runs | awk '{FS=" \n\t"}{printf "%d", $1}' | head -n 1` 
    if test x$run_for_acc = "x" || test $run_for_acc -lt 1; then 
	echo "No run for acceptance correction specified" > /dev/stderr 
	exit 1
    fi

   now=`date '+%Y%m%d_%H%M'` 
   outputs

    # Write settings to a file, which we later can source 
    cat > ${dotconf} <<EOF
name="$name"
runs=${runs}
sys=$sys
snn=$snn
field=$field
real_dir=${real_dir}
real_pat=${real_pat}
mc_dir=${mc_dir}
mc_pat=${mc_pat}
my_real_dir=${my_real_dir}
my_mc_dir=${my_mc_dir}
par=${par}
now=${now}
EOF

    if test $noact -lt 1 ; then 
	mkdir -p ${name}_acc_${now}
    fi
    echo "Make acceptance corrections" 
    (cd ${name}_acc_${now} && \
	accGen $run_for_acc && \
	upload )
}    

# --- Run set-ups ----------------------------------------------------
cleanup()
{
    rm -rf \
	${name}_acc_${now} \
	${name}_mccorr_${now} \
	${name}_mceloss_${now} \
	${name}_eloss_${now} \
	${name}_mcaod_${now} \
	${name}_aod_${now} \
	${name}_mcdndeta_${now} \
	${name}_dndeta_${now} \
	${name}_corrs_${now}
}

# --- Check settings -------------------------------------------------
check()
{
   if test "x$runs" = "x" || test ! -f $runs ; then 
       echo "List of run file $runs not found" > /dev/stderr 
       exit 1
   fi
   if test "X$name" = X ; then 
       echo "No name specified" > /dev/stderr 
       exit 1
   fi
   if test "x$sys" = "x" ; then 
       echo "No collision system specified" > /dev/stderr 
       exit 1
   fi
   if test "x$snn" = "x" ; then 
       echo "No center of mass energy specified" > /dev/stderr 
       exit 1
   fi
   if test "x$field" = "x" ; then 
       echo "No L3 field setting specified" > /dev/stderr 
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
   if test "X$1" != "Xsetup" && test x$now = x ; then 
       echo "No date/time specified" > /dev/stderr 
       exit 1
   fi
   case $sys in 
       pp|p-p)            sys=1 ;; 
       pbpb|pb-pb|aa|a-a) sys=2 ;; 
       ppb|p-pb|pa|p-a)   sys=3 ;;
       1|2|3)                   ;; 
       *) echo "$0: Unknown system: $sys" ; exit 1 ;;
   esac

   cat <<EOF
Name:			$name
Run file:		${runs}
Collision system:	$sys
sqrt(s_NN):		${snn}GeV
L3 Field:		${field}kG
Real input directory:	${real_dir}
Real file pattern:	${real_pat}
MC input directory:	${mc_dir}
MC file pattern:	${mc_pat}
Real output:		${my_real_dir}
MC output directory:	${my_mc_dir}
Use PAR files:		${par}
Date & time:            ${now}

EOF
}

# --- Run the train --------------------------------------------------
# Usage:
# 
allAboard()
{
    type=$1 ; shift 
    cl=
    nme=${name}_${type}
    tree=esdTree
    opts=""
    uopt=""
    mc=0
    dir=$real_dir
    pat=$real_pat

    case $type in 
	mc*) mc=1 ;; 
	*) ;;
    esac
    case $type in 
	*corr)  cl=MakeMCCorrTrain ; mc=1 ;;
	*eloss) cl=MakeFMDELossTrain ;;  
	*aod)   cl=MakeAODTrain 
	    opts="--corr=../${name}_corrs_${now} --sys=${sys} --snn=${snn} --field=${field}"
	    ;;
	*dndeta) cl=MakedNdetaTrain 
	    tree=aodTree 
	    uopt="&concat"
	    opts="${opts}"
	    ;;
	*) echo "$0: Unknown type of train: $type" > /dev/stderr ; exit 1 ;;
    esac
    case $type in 
	*corr|*eloss|*aod)
	    if test $mc -gt 0; then 
		uopt="${uopt}&mc"
		dir=$mc_dir
		pat=$mc_pat
	    fi
	    ;;
	*dndeta)
	    dir=$my_real_dir
	    pat="*/AliAOD.root"
	    if test $mc -gt 0 ; then 
		uopt="${uopt}&mc"
		dir=$my_mc_dir
	    fi
	    ;;
    esac
    case $type in 
	*aod|*dndeta) 
	    if test $sys -gt 1 ; then opts="${opts} --cent" ; fi ;; 
	*)
	    ;;
    esac
    if test $par -gt 0 ; then 
	uopt="${uopt}&par"
    fi
    url="alien://${dir}?run=${runs}&pattern=${pat}${uopt}${aliroot}${root}#${tree}"
    opts="${opts} --date=${now} --class=$cl --name=$nme"
    
    echo "Running train: runTrain2 ${opts} --url=${url} $@" 
    if test $noact -gt 0 ; then return ; fi

    runTrain2 ${opts} --overwrite --url=${url} $@ 

    cat <<EOF
Check https://alimonitor.cern.ch/users/jobs.jsp for progress 

Remember to do 

  (cd ${nme}_${now} && aliroot -l -b -q Terminate.C)

until the final merge stage, and then do 

  (cd ${nme}_${now} && aliroot -l -b -q Download.C) 

to get the results. 
EOF
    case $type in 
	*corr|*esd) 
	    cat <<EOF
Then, do 

  (cd ${nme}_${now} && aliroot -l -b -q Extract.C)
  (cd ${nme}_${now} && aliroot -l -b -q 'Upload.C("local://${here}/${name}_corrs_${now}")')

to upload the results to our local corrections store. 
EOF
	    ;; 
	*aod)
	    cat <<EOF
Then, do 
 
  (cd ${nme}_${now} && aliroot -l ${fwd_dir}/DrawAODSummary.C)

to get a PDF of the diagnostics histograms
EOF
	    ;;
	*dndeta)
	    cat <<EOF
Then, do 
 
  (cd ${nme}_${now} && aliroot -l ${fwd_dir}/DrawdNdetaSummary.C)

to get a PDF of the diagnostics histograms, and 

  (cd ${nme}_${now} && aliroot -l draw.C)

to get the final plot. 
EOF
	    ;;
    esac
}


# === Wrappers =======================================================
# --- Run all correction jobs ----------------------------------------
corrs()
{
    allAboard mccorr $@
    allAboard mceloss $@
    allAboard eloss $@
}
corrs_terminate() 
{
    (cd ${name}_mccorr_${now}  && terminate)
    (cd ${name}_mceloss_${now} && terminate)
    (cd ${name}_eloss_${now}   && terminate)
}
corrs_upload() 
{
    (cd ${name}_mccorr_${now}  && extract_upload)
    (cd ${name}_mceloss_${now} && extract_upload)
    (cd ${name}_eloss_${now}   && extract_upload)
}
corrs_draw()
{
    echo "Draw does not make sense for Correction jobs"
}
# --- Run all AOD jobs -----------------------------------------------
aods()
{
    allAboard mcaod $@
    allAboard aod $@
}
aods_terminate() 
{
    (cd ${name}_mcaod_${now} && terminate)
    (cd ${name}_aod_${now}   && terminate)
}
aods_upload()
{
    echo "Upload does not make sense for AOD jobs"
}
aods_draw() 
{
    (cd ${name}_mcaod_${now} && draw ${fwd_dir}/DrawAODSummary.C)
    (cd ${name}_aod_${now}   && draw ${fwd_dir}/DrawAODSummary.C)
}

# --- Run all dN/deta jobs -------------------------------------------
dndetas()
{
    allAboard mcdndeta $@
    allAboard dndeta $@
}
dndetas_terminate() 
{
    (cd ${name}_mcdndeta_${now} && terminate)
    (cd ${name}_dndeta_${now}   && terminate)
}
dndetas_upload()
{
    echo "Upload does not make sense for dN/deta jobs"
}
dndetas_draw() 
{
    (cd ${name}_mcdndeta_${now} && draw ${fwd_dir}/DrawdNdetaSummary.C && \
	script draw.C)
    (cd ${name}_dndeta_${now}   && draw ${fwd_dir}/DrawdNdetaSummary.C && \
	script draw.C)
}

# === Executable code
# --- Source settings if found ---------------------------------------
if test -f $dotconf ; then 
    source $dotconf 
fi


# --- Process command line -------------------------------------------
what=
step=
while test $# -gt 0 ; do
    arg=$1 
    opt=
    case $1 in 
	--) shift ; break ;;
	--*=*) 
	    arg=`echo $1 | sed 's/=.*//'` ;
	    opt=`echo $1 | sed 's/--[^=][^=]*=//'` 
	    ;;
	--*)
	    ;;
	-h|-N|-H|-a) ;;
	-*) opt=$2 ; shift ;; 
    esac
    shift 

    case $arg in 
	-r|--runs)         runs=$opt     ;; 
	-n|--name)         name=$opt     ;; 
	-S|--sys)          sys=`echo $opt | tr '[A-Z]' '[a-z]'` ;;
	-E|--snn)          snn=$opt      ;; 
	-F|--field)        field=$opt    ;;
	-d|--real-dir)     real_dir=$opt ;;
	-p|--real-pattern) real_pat=$opt ;; 
	-D|--mc-dir)       mc_dir=$opt   ;; 
	-P|--mc-pattern)   mc_pat=$opt   ;; 
	-w|--what)         what=`echo $opt | tr '[A-Z]' '[a-z]'`  ;; 
	-s|--step)         step=`echo $opt | tr '[A-Z]' '[a-z]'`  ;; 
	-N|--noact)        noact=1       ;;
	-a|--par)          par=1         ;;
	-h|--help)         usage         ; exit 0 ;; 
	-H|--manual)       manual        ; exit 0 ;;
	*) echo "$0: Unknown option $arg"  ; exit 1 ;; 
    esac
done 

# --- Check settings -------------------------------------------------
check $what

# --- Select what to do ----------------------------------------------
func=
case $what in 
    setup)   setup ; exit 0 ;; 
    clean)   cleanup ; exit 0 ;;
    corr*)   func=corrs;; 
    aod*)    func=aods ;; 
    dndeta*) func=dndeta ;; 
    *) echo "$0: Unknown operation: $what" > /dev/stderr ; exit 1 ;;
esac
    
case x$step in 
    x|xfull) ;; 
    xterm*) func=${func}_terminate ;; 
    xup*)   func=${func}_upload ;; 
    xdr*)   func=${func}_draw ;;
    *) echo "$0: Unknown step $step" > /dev/stderr ; exit 1 ;;
esac

echo "Will execute $func" 
$func $@ 

#
# EOF
#

