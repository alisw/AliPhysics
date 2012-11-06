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
#   $0 --what setup  \
#       --name=LHC10h --sys=pbpb --snn=2760 --field=5 \
#       --real-dir=/alice/data/2010/LHC10h \
#       --real-pattern=ESDs/pass2/*/AliESDs.root \
#       --mc-dir=/alice/sim/LHC11a7 \
#       --mc-pattern=*/AliESDs.root \
#       --runs=runs.list --par
# 
# Note, all the settings are written to the file .config in the
# current directory, so you do not need to give the parameters at
# subsequent steps.  As an alternative to giving the parameters, one
# can create the file by hand.
# 
# Next, we need to generate the corrections.  Do 
# 
#   $0 --what corrs 
# 
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what corrs --step terminate 
# 
# enough times to get the final merged result.  Next, we need to download 
# the results
# 
#   $0 --what corrs --step download 
# 
# and extract and upload the corrections to our local corrections folder
# 
#   $0 --what corrs --step upload 
# 
# Now we can submit our AOD generation jobs.  Do 
# 
#   $0 --what aods 
# 
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what aods --step terminate 
# 
# enough times to get the final merged result.  Next, we need to download 
# the results
# 
#   $0 --what aods --step download 
# 
# and we can draw the summary results 
# 
#   $0 --what aods --step draw 
# 
# Now, we should do the dN/deta analysis.  Do 
# 
#   $0 --what dndetas
# 
# and then monitor output on MonAlisa.  When enough has finished, execute 
# 
#   $0 --what dndetas --step terminate 
# 
# enough times to get the final merged result.  Next, we need to download 
# the results
# 
#   $0 --what dndetas --step download 
# 
# and we can draw the final results 
# 
#   $0 --what dndetas --step draw 
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

real_dir=
real_pat=
mc_dir=
mc_pat=
my_real_dir=/alice/cern.ch/users/
my_mc_dir=

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
  download    Download results 
  upload      Upload corrections (only for OPERATION=corrs)
  draw        Draw (partial) results
EOF
}

# --- Manual ---------------------------------------------------------
manual()
{
    grep ^# $0 | \
	sed -n -e '/BEGIN_MANUAL/,/END_MANUAL/ p' | \
	sed -e '/\(BEGIN\|END\)_MANUAL/ d' -e 's/^# //' -e "s,\$0,$0,"
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
   
# --- Run script -----------------------------------------------------
script()
{
    scr=$1 ; shift 
    args=$1 
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
    script $ALICE_ROOT/PWGLF/FORWARD/analysis2/corrs/ExtractAcceptance.C \
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
    script Download.C 
}
# --- Extract corrections --------------------------------------------
extract()
{
    script Extract.C 
}
# --- Upload a file --------------------------------------------------
upload()
{
    script Upload.C \"file://${here}/${corrs}\"
}
# --- Extract and upload ---------------------------------------------
extract_upload()
{
  for i in *.zip ; do 
      if test ! -f $i ; then continue ; fi 
      echo "Extracting and uploading from $i"
      if test $noact -gt 0 ; then continue ; fi
      rm -rf tmp
      (mkdir -p tmp && \
	  cd tmp && \ 
	  unzip ../$i && \
	      script ../Extract.C && 
	  upload)
  done
}

# --- Run set-ups ----------------------------------------------------
cleanup()
{
    rm -rf ${name}_acc \
	${name}_mccorr \
	${name}_mceloss \
	${name}_eloss \
	${name}_mcaod \
	${name}_aod \
	${name}_mcdndeta \
	${name}_dndeta \
	tmp \
	${corrs}
}

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
	mkdir -p ${name}_acc
    fi
    echo "Make acceptance corrections" 
    (cd ${name}_acc && \
	accGen $run_for_acc && \
	upload )
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
	    opts="--corr=../${corrs} --sys=${sys} --snn=${snn} --field=${field}"
	    ;;
	*dndeta) cl=MakedNdetaTrain 
	    tree=aodTree 
	    uopt="&concat"
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
  (cd ${nme}_${now} && aliroot -l -b -q 'Upload.C("local://${here}/${corrs}")')

to upload the results to our local corrections store. 
EOF
	    ;; 
	*aod)
	    cat <<EOF
Then, do 
 
  (cd ${nme}_${now} && aliroot -l $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawAODSummary.C)

to get a PDF of the diagnostics histograms
EOF
	    ;;
	*dndeta)
	    cat <<EOF
Then, do 
 
  (cd ${nme}_${now} && aliroot -l $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawdNdetaSummary.C)

to get a PDF of the diagnostics histograms, and 

  (cd ${nme}_${now} && aliroot -l draw.C)

to get the final plot. 
EOF
	    ;;
    esac
}


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
corrs_download() 
{
    (cd ${name}_mccorr_${now}  && download)
    (cd ${name}_mceloss_${now} && download)
    (cd ${name}_eloss_${now}   && download)
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
aods_download() 
{
    (cd ${name}_mcaod_${now} && download)
    (cd ${name}_aod_${now}   && download)
}
aods_upload()
{
    echo "Upload does not make sense for AOD jobs"
}
aods_draw() 
{
    (cd ${name}_mcaod_${now} && \
	script $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawAODSummary.C)
    (cd ${name}_aod_${now} && \
	script $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawAODSummary.C)
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
dndetas_download() 
{
    (cd ${name}_mcdndeta_${now} && download)
    (cd ${name}_dndeta_${now}   && download)
}
dndetas_upload()
{
    echo "Upload does not make sense for dN/deta jobs"
}
dndetas_draw() 
{
    (cd ${name}_mcdndeta_${now} && \
	script $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawdNdetaSummary.C && \
	script draw.C)
    (cd ${name}_dndeta_${now} && \
	script $ALICE_ROOT/PWGLF/FORWARD/analysis2/DrawdNdetaSummary.C && \
	script draw.C)
}

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
    
case $step in 
    full) ;; 
    term*) func=${func}_terminate ;; 
    down*) func=${func}_download ;; 
    up*)   func=${func}_upload ;; 
    dr*)   func=${func}_draw ;;
esac

echo "Will execute $func" 
$func $@ 

#
# EOF
#

