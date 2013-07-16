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
# enough times to get the final merged result.  Next, we need to download 
# the results and we can draw the summary and final plot
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

runs=
mcruns=
name=
sys=0
snn=0
field=0
corrs=
dotconf=.config
here=${PWD}
par=0
noact=0
# aliroot="&aliroot=v5-03-75pATF-AN"
# root="root=v5-34-02-1"
fwd_dir=$ALICE_ROOT/PWGLF/FORWARD/analysis2

real_dir=
real_pat=
mc_dir=
mc_pat=
my_real_dir=
my_mc_dir=
uuopts=
watch=0

# === Various functions ==============================================
# --- Usage ----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 --what OPERATION [OPTIONS]

Options:
  -r,--runs=FILENAME        Specify list of runs file ($runs)
  -R,--mc-runs=FILENAME     Specify list of MC runs file ($mcruns)
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
  -c,--corrections=DIR      Directory where corrections are stored ($corrs)
  -W,--watch                Watch for job status and terminate automatically
  -a,--par                  Use par files ($par)
  -u,--url-opts=OPTIONS	    Additional user options ($uuopts)
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
    prog=`basename $0`
    grep ^# $0 | \
	sed -n -e '/BEGIN_MANUAL/,/END_MANUAL/ p' | \
	sed -e '/\(BEGIN\|END\)_MANUAL/ d' -e 's/^# //' \
	    -e "s,\$0,$prog,"
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
    script ${fwd_dir}/corrs/ExtractAcceptance.C "${run}"
}

# --- Extract corrections --------------------------------------------
terminate()
{
    script Terminate.C 
}
# --- Extract corrections --------------------------------------------
download()
{
    if test -f .download ; then 
	echo "Already downloaded in `basename $PWD`"
	return 0 
    fi
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
# --- Extract corrections --------------------------------------------
extract_up()
{
    if test -f .extract ; then 
	echo "Aldready extracted in `basename $PWD`" 
	return 0;
    fi 
    echo "= Extracting"
    script ../Extract.C > /dev/null 2>&1
    touch .extract
}
# --- Upload a file --------------------------------------------------
upload()
{
    if test -f .upload ; then 
	echo "Already uploaded in `basename $PWD`"
	return 0 
    fi
    echo "= Uploading"
    script Upload.C \"file://${here}/${name}_corrs_${now}/\" >/dev/null 2>&1
    touch .upload 
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
	extract_up 
	upload
	cd ..
    done
}

# --- Draw -----------------------------------------------------------
draw()
{
    local scr=$1 
    download 
    for i in *.zip ; do 
	if test "X$i" = "X*.zip" ; then continue ; fi
	echo "Will extract $i"
	d=`basename $i .zip` 
	if test ! -d $d ; then 
	    mkdir -p $d 
	    unzip $i -d $d
	fi
	(cd $d && script $scr)
    done
}
dndeta_draw()
{
    local d=$1 
    (cd $d && \
	draw ${fwd_dir}/DrawdNdetaSummary.C && \
	draw ../Draw.C)
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

# === Trains =========================================================
# --- Run set-ups ----------------------------------------------------
setup()
{
    run_for_acc=`grep -v ^# $runs | awk '{FS=" \n\t"}{printf "%d\n", $1}' | head -n 1` 
    if test x$run_for_acc = "x" || test $run_for_acc -lt 1; then 
	echo "No run for acceptance correction specified" > /dev/stderr 
	exit 1
    fi

   now=`date '+%Y%m%d_%H%M'` 
   outputs

   # Write settings to a file, which we later can source 
   dumpvar=
   if test $par -gt 0 ; then dumpvar="--par " ; fi 
   if test $watch -gt 0 ; then dumpvar="${dumpvar} --watch " ; fi 
   cat > ${dotconf} <<EOF
# Settings:
name="$name"
runs=${runs}
mcruns=${mcruns}
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
watch=${watch}
uuopts="${uuopts}"
# Options
if false ; then 
  $0 --what=setup --name="$name" --runs="$runs" --mcruns="$mcruns" \
  --sys="$sys" --snn="$snn" --field="$field" \
  --real-dir="${real_dir}" --real-pattern="${real_pat}" \
  --mc-dir="${mc_dir}" --mc-pattern="${mc_pat}" \
  --now=${now} --url-opts="${uuopts}" ${dumpvar}
fi
EOF
   corrdir=${name}_corrs_${now}
   if test "x$corrs" != "x" && test -d ${corrs} ; then 
       echo "Linking ${corrs} to ${corrdir}"
       ln -sf $corrs ${corrdir}
       ln -sf $corrs last_${name}_corrs	       
       corrdir=$corrs
   elif test $noact -lt 1 ; then 
	mkdir -p ${name}_acc_${now}
	mkdir -p ${corrdir}
	rm -f last_${name}_acc last_${name}_corrs
	ln -sf ${name}_acc_${now} last_${name}_acc
	ln -sf ${name}_corrs_${now} last_${name}_corrs
	cat <<-EOF > ${corrdir}/Browse.C
 	TObject* Browse()
	{
	  const char* fwd = "$ALICE_ROOT/PWGLF/FORWARD/analysis2";
	  if (!gROOT->GetClass("AliOADBForward"))
	    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
	  gROOT->LoadMacro(Form("%s/corrs/ForwardOADBGui.C++g", fwd));
	  
	  AliOADBForward* db = new AliOADBForward;
	  db->Open("fmd_corrections.root", "*");
	  
	  ForwardOADBGui(db);

	  return db;
	}
	EOF
	echo "Make acceptance corrections" 
	(cd ${name}_acc_${now} && \
	    accGen $run_for_acc && \
	    upload )
   fi
   for i in fmd_corrections.root spd_corrections.root deadstrips.C ; do 
       if test ! -f ${corrdir}/$i ; then continue ; fi 
       echo "Linking ${corrdir}/$i here"
       ln -fs ${corrdir}/$i . 
   done
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
	${name}_mcdndeta_*${now} \
	${name}_dndeta_*${now} \
	${name}_corrs_${now} \
	last_${name}_acc \
	last_${name}_corrs \
	last_${name}_mceloss \
	last_${name}_eloss \
	last_${name}_mcaod \
	last_${name}_aod \
	last_${name}_mcdndeta* \
	last_${name}_dndeta* \
	build.log 
    if test -L fmd_corrections.root ; then 
	rm fmd_corrections.root
    fi
    if test -L spd_corrections.root ; then 
	rm spd_corrections.root
    fi
}

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

# --- Check settings -------------------------------------------------
check()
{
    local w=$1
    if test "x$runs" = "x" || test ! -f $runs ; then 
	echo "List of run file $runs not found" > /dev/stderr 
	exit 1
    fi
    if test "x$mcruns" = "x" ; then mcruns=$runs ; fi 
    if test ! -f $mcruns ; then 
	echo "List of MC runs file $mcruns not found" > /dev/stderr 
	exit 1
    fi
    if test "X$name" = X ; then 
	echo "No name specified" > /dev/stderr 
	exit 1
    fi
    # if test "x$sys" = "x" ; then 
    #     echo "No collision system specified" > /dev/stderr 
    #     exit 1
    # fi
    # if test "x$snn" = "x" ; then 
    #     echo "No center of mass energy specified" > /dev/stderr 
    #     exit 1
    # fi
    # if test "x$field" = "x" ; then 
    #     echo "No L3 field setting specified" > /dev/stderr 
    #     exit 1
    # fi
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
    if test "X$w" != "Xsetup" && test "x$now" = "x" ; then 
	echo "No date/time specified" > /dev/stderr 
	exit 1
    fi
    # sys==0 is OK - autoselect
    case x$sys in 
        xpp|xp-p)              sys=1 ;; 
	xpbpb|xpb-pb|xaa|xa-a) sys=2 ;; 
	xppb|xp-pb|xpa|xp-a)   sys=3 ;;
	x0|x1|x2|x3)                 ;; 
	x)                     sys=0 ;;
	*) echo "$0: Unknown system: $sys" ; exit 1 ;;
    esac

    check_token
}

# --- Show the setup -------------------------------------------------
print_setup()
{
    cat <<EOF
Name:			$name
Run file:		${runs}
MC Run file:            ${mcruns}
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
Additional URL options: ${uuopts}
EOF
}

# --- Run the train --------------------------------------------------
# Usage:
# 
allAboard()
{
    type=$1 ; shift 
    trig=$1 ; shift
    cl=
    nme=${name}_${type}
    tree=esdTree
    opts="--batch"
    uopt="&merge=50&split=50&aliroot=last,regular"
    mc=0
    dir=$real_dir
    pat=$real_pat
    rl=$runs

    case $type in 
	mc*) 
	    mc=1  
	    # Default dirs are production dirs 
	    dir=$mc_dir
	    pat=$mc_pat
	    rl=$mcruns
	    ;;
	*) ;;
    esac
    case $type in 
	*corr)  cl=MakeMCCorrTrain ; mc=1 ;;
	*eloss) cl=MakeFMDELossTrain ;;  
	*aod)   cl=MakeAODTrain 
	    opts="${opts} --corr=."
	    # opts="--corr=${name}_corrs_${now} --cent"
	    # if test $sys -gt 0 && test $snn -gt 0 ; then 
	    # 	opts="$opts --sys=${sys} --snn=${snn} --field=${field}"
	    # fi
	    ;;
	*dndeta) cl=MakedNdetaTrain 
	    tree=aodTree 
	    uopt="${uopt}&concat"
	    opts="${opts} --cut-edges"
	    case x$trig in 
		xinel)    
		    opts="$opts --scheme=trigger,event,background}" 
		    opts="$opts --trig=INEL" 
		    ;;
		xnsd)     
		    opts="$opts --scheme=trigger,event"
		    opts="$opts --trig=V0AND"
		    ;;
		xinelgt0) 
		    opts="$opts --scheme=trigger,event"
		    opts="$opts --trig=INELGT0"
		    ;;
		x*) trig= ;;
	    esac
	    if test "x$trig" != "x" ; then 
		nme="${nme}_${trig}"
	    fi
	    # Modify for input dir for our files
	    dir=$my_real_dir
	    pat="*/AliAOD.root"
	    if test $mc -gt 0 ; then 
		dir=$my_mc_dir
		opts="$opts --mc"
	    fi
	    ;;
	*) echo "$0: Unknown type of train: $type" > /dev/stderr ; exit 1 ;;
    esac
    # add centrality flag if we do not know what collision system we're 
    # looking at, or it's PbPb or pPb. 
    case $sys in 
	0|2|3) opts="$opts --cent" ;; 
	1)                         ;;
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
    url="alien://${dir}?run=${rl}&pattern=${pat}${uopt}${aliroot}${root}#${tree}"
    opts="${opts} --include=$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains"
    opts="${opts} --date=${now} --class=$cl --name=$nme --verbose=0"
    
    echo "Running train: runTrain ${opts} --url=${url} $@" 
    if test $noact -gt 0 ; then return ; fi

    runTrain ${opts} --overwrite --url=${url} $@ 

    if test $watch -lt 1 ; then 
	cat <<-EOF
	Check https://alimonitor.cern.ch/users/jobs.jsp for progress 
	
	Remember to do 
	
	  (cd ${nme}_${now} && aliroot -l -b -q Terminate.C)
	
	until the final merge stage, and then do 
		
	  (cd ${nme}_${now} && aliroot -l -b -q Download.C) 
	
	to get the results. 
	EOF
	
	case $type in 
	    *corr|*esd) 
		cat <<-EOF
	Then, do 
	
	  (cd ${nme}_${now} && aliroot -l -b -q Extract.C)
	  (cd ${nme}_${now} && aliroot -l -b -q 'Upload.C("local://${here}/${name}_corrs_${now}")')
	
	to upload the results to our local corrections store. 
	EOF
		;; 
	    *aod)
		cat <<-EOF
	Then, do 
	 
	  (cd ${nme}_${now} && aliroot -l ${fwd_dir}/DrawAODSummary.C)
	
	to get a PDF of the diagnostics histograms
	EOF
		;;
	    *dndeta)
		cat <<-EOF
	Then, do 
	 
	  (cd ${nme}_${now} && aliroot -l ${fwd_dir}/DrawdNdetaSummary.C)
	
	to get a PDF of the diagnostics histograms, and 
	
	  (cd ${nme}_${now} && aliroot -l Draw.C)
	
	to get the final plot. 
	EOF
		;;
	esac
    else 
	echo "Now waiting for jobs to finish"
	(cd ${nme}_${now} && \
	    nice aliroot -l -b -x -q Watch.C\(1\) > watch.log 2>&1 &)
    fi
}


# === Wrappers =======================================================
# --- Run all correction jobs ----------------------------------------
corrs()
{
    allAboard mccorr "" $@
    allAboard mceloss "" $@
    allAboard eloss "" $@
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
    rm -f fmd_corrections.root spd_corrections.root 
    ln -s ${name}_corrs_${now}/fmd_corrections.root .
    ln -s ${name}_corrs_${now}/spd_corrections.root .
}
corrs_draw()
{
    echo "Draw does not make sense for Correction jobs"
}
# --- Run all AOD jobs -----------------------------------------------
aods()
{
    allAboard mcaod "" $@
    allAboard aod "" $@
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
    if test $sys -eq 1 ; then 
	allAboard dndeta inel $@
	allAboard dndeta nsd $@
	allAboard dndeta inelgt0 $@
    else
	allAboard dndeta "" $@
    fi
}
dndetas_terminate() 
{
    (cd ${name}_mcdndeta_${now} && terminate)
    if test $sys -eq 1 ; then 
	(cd ${name}_dndeta_inel_${now}      && terminate)
	(cd ${name}_dndeta_nsd_${now}       && terminate)
	(cd ${name}_dndeta_inelgt0_${now}   && terminate)
    else
	(cd ${name}_dndeta_${now}   && terminate)
    fi
}
dndetas_upload()
{
    echo "Upload does not make sense for dN/deta jobs"
}
dndetas_draw() 
{
    dndeta_draw ${name}_mcdndeta_${now}
    if test $sys -eq 1 ; then 
	dndeta_draw ${name}_dndeta_inel_${now} 
	dndeta_draw ${name}_dndeta_nsd_${now} 
	dndeta_draw ${name}_dndeta_inelgt0_${now} 
    else
	dndeta_draw ${name}_dndeta_${now} 
    fi
}

# === Procedual code =================================================
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
	-R|--mc-runs)      mcruns=$opt   ;;
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
	-W|--watch)        let watch=\!$watch ;;
	-N|--noact)        noact=1            ;;
	-c|--corrections)  corrs=$opt    ;;
	-a|--par)          par=1         ;;
	-u|--url-opts)     uuopts="$opt" ;;
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
    dndeta*) func=dndetas ;; 
    *) echo "$0: Unknown operation: $what" > /dev/stderr ; exit 1 ;;
esac
print_setup

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

