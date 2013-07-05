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
# options to the train can be passed after --.  if you passed the
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

run=
name=
sys=0
snn=0
field=0
corrs=
dotconf=.config
here=${PWD}
par=0
noact=0
nwrks=0
fwd_dir=$ALICE_ROOT/PWGLF/FORWARD/analysis2

real_dir=
real_pat=
real_idx=
mc_dir=
mc_pat=
mc_idx=
my_real_dir=
my_mc_dir=
uuopts=
# Latest: 0.91
inel_eff=1
# Latest: 0.94
nsd_eff=1
inelgt0_eff=1

# === Various functions ==============================================
# --- Usage ----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 --what OPERATION [OPTIONS]

Options:
  -r,--run=NUMBER           Specify run number ($run)
  -n,--name=STRING          Base name of jobs ($name)
  -S,--sys=SYSTEM           Collision system ($sys)
  -E,--snn=ENERGY           Center of mass energy per nuclean pair ($snn)
  -F,--field=FIELD          L3 magnetic field ($field)
  -d,--real-dir=DIR         Directory holding real data ($real_dir)
  -p,--real-pattern=PATTERN Glob pattern to match when searching ($real_pat)
  -D,--mc-dir=DIR           Directory holding MC data ($mc_dir)
  -P,--mc-pattern=PATTERN   Glob pattern to match when searching ($mc_pat)
  -s,--step=STEP            Run stage ($step)
  -w,--what=TRAINS          What to do 
  -c,--corrections=DIR      Directory where corrections are stored ($corrs)
  -W,--workers=N            Number of workers ($nwrks)
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
  upload      Upload corrections (only for TRAINS=corrs)
  draw        Draw (partial) results
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
    check_token
    local run=$1 
    script ${fwd_dir}/corrs/ExtractAcceptance.C "${run}"
}
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
    if test ! -d $d || test ! -f $d ;then 
	echo "Specified input for $t data is not a directory or file " \
	    > /dev/stderr 
	return;
    fi
    if test $m -gt 0 ; then 
	n="&mc"
    fi
    aliroot -l -b <<EOF
.L $ALICE_ROOT/PWGLF/FORWARD/trains/ChainBuilder.C++
ChainBuilder::CreateCollection("${o}", "file://${d}?recursive&scan&pattern=${p}${n}#esdTree");
.q
EOF
    
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
    extract
    upload
}

# --- Draw -----------------------------------------------------------
draw()
{
    script $@
}
dndeta_draw()
{
    local d=$1 
    echo "=== $d ================================================"
    (cd $d && \
	draw ${fwd_dir}/DrawdNdetaSummary.C && \
	draw Draw.C)
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

# === Trains =========================================================
# --- Run set-ups ----------------------------------------------------
setup()
{
    if test x$run = "x" || test $run -lt 1; then 
	echo "No run for acceptance correction specified" > /dev/stderr 
	exit 1
    fi

   now=`date '+%Y%m%d_%H%M'` 
   outputs

   # Write settings to a file, which we later can source 
   dumpvar=
   if test $par -gt 0 ; then dumpvar="--par " ; fi 
   cat > ${dotconf} <<EOF
# Settings:
name="$name"
run=${run}
sys=$sys
snn=$snn
field=$field
real_dir=${real_dir}
real_pat=${real_pat}
real_idx=${real_idx}
mc_dir=${mc_dir}
mc_pat=${mc_pat}
mc_idx=${mc_idx}
my_real_dir=${my_real_dir}
my_mc_dir=${my_mc_dir}
par=${par}
now=${now}
uuopts="${uuopts}"
# Trigger efficiencies - edit here to set them 
inel_eff=$inel_eff
inelgt0_eff=$inelgt0_eff
nsd_eff=$nsd_eff
# Options
if false ; then 
  $0 --what=setup --name="$name" --run="$run" \
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
	${name}_corrs_${now} \
	${real_idx} \
	${mc_idx} 
}


# --- Check settings -------------------------------------------------
check()
{
    local w=$1
    if test "x$run" = "x" || test $run -lt 1 ; then 
	echo "Run not specified, or invalid ($run)" > /dev/stderr 
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
    
    ncpu=`cat /proc/cpuinfo|sed -n 's/^processor[ \t]*: \(.*\)/\1/p'|wc -l`
    if test "x$nwrks" = "x" || test $nwrks -lt 1 ; then 
	let nwrks=7*$ncpu/10
	echo "Setting number of workers to $nwrks / $ncpu"
    fi
}

# --- Show the setup -------------------------------------------------
print_setup()
{
    cat <<EOF
Name:			${name}
Run:			${run}
Collision system:	${sys}
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
Number of workers:      ${nwrks}/${ncpu}
Trigger efficiencies:   
  INEL:                 ${inel_eff}
  INEL>0:               ${inelgt0_eff}
  NSD:                  ${nsd_eff}
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
    opts=""
    # opts="--batch"
    uopt="mode=default&workers=${nwrks}"
    mc=0
    inp=${real_idx}
    # dir=$real_dir
    # pat=$real_pat
    rl=$runs

    case $type in 
	mc*) 
	    mc=1  
	    # Default dirs are production dirs 
	    inp=${mc_idx}
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
	    opts="${opts} --cut-edges"
	    case x$trig in 
		xinel)    
		    opts="$opts --scheme=trigger,event,background" 
		    opts="$opts --trig=INEL --trigEff=$inel_eff" 
		    ;;
		xnsd)     
		    opts="$opts --scheme=trigger,event"
		    opts="$opts --trig=V0AND --trigEff=$nsd_eff"
		    ;;
		xinelgt0) 
		    opts="$opts --scheme=trigger,event"
		    opts="$opts --trig=INELGT0 --trigEff=$inelgt0_eff"
		    ;;
		x*) trig= ;;
	    esac
	    if test "x$trig" != "x" ; then 
		nme="${nme}_${trig}"
	    fi
	    # Modify for input dir for our files
	    inp=$my_real_dir/AliAOD.root
	    if test $mc -gt 0 ; then 
		inp=$my_mc_dir/AliAOD.root
		opts="$opts --mc"
	    fi
	    ;;
	*multdists) 
	    cl=MakeMultDistsTrain 
	    tree=aodTree 
	    # Modify for input dir for our files
	    inp=$my_real_dir/AliAOD.root
	    if test $mc -gt 0 ; then 
		inp=$my_mc_dir/AliAOD.root
	    fi
	    ;;
	*) echo "$0: Unknown type of train: $type" > /dev/stderr ; exit 1 ;;
    esac
    if test ! -f $inp ; then 
	echo "No input for $nme, giving up" > /dev/stderr 
	return 
    fi
    # add centrality flag if we do not know what collision system we're 
    # looking at, or it's PbPb or pPb. 
    case $sys in 
	1)                     ;;
	*) opts="$opts --cent" ;; 
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
    opts="${opts} --include=$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains"
    opts="${opts} --date=${now} --class=$cl --name=$nme --verbose=0"
    
    echo "Running train: runTrain ${opts} --url=${url} $@" 
    if test $noact -gt 0 ; then return ; fi

    runTrain ${opts} --url=${url} $@ 
}


# === Wrappers =======================================================
# --- Run all correction jobs ----------------------------------------
corrs()
{
    allAboard mccorr "" $@
    allAboard mceloss "" $@
    allAboard eloss "" $@
}
corrs_upload() 
{
    if test "X$mc_dir" != "X" ; then 
	(cd ${name}_mccorr_${now}  && extract_upload)
	(cd ${name}_mceloss_${now} && extract_upload)
    fi
    (cd ${name}_eloss_${now}   && extract_upload)
    rm -f fmd_corrections.root spd_corrections.root 
    ln -s ${name}_corrs_${now}/fmd_corrections.root .
    ln -s ${name}_corrs_${now}/spd_corrections.root .
}
corrs_draw()
{
    if test "X$mc_dir" != "X" ; then 
	(cd ${name}_mccorr_${now}  && draw ${fwd_dir}/DrawMCCorrSummary.C)
	(cd ${name}_mceloss_${now} && draw ${fwd_dir}/corrs/DrawCorrELoss.C 1)
    fi
    (cd ${name}_eloss_${now}   && draw ${fwd_dir}/corrs/DrawCorrELoss.C 0)
}
# --- Run all AOD jobs -----------------------------------------------
aods()
{
    allAboard mcaod "" $@
    allAboard aod "" $@
}
aods_upload()
{
    echo "Upload does not make sense for AOD jobs"
}
aods_draw() 
{
    (cd ${name}_mcaod_${now} && draw Summarize.C)
    (cd ${name}_aod_${now}   && draw Summarize.C)
}

# --- Run all dN/deta jobs -------------------------------------------
dndetas()
{
    if test $sys -eq 1 ; then 
	allAboard mcdndeta inel    $@
	allAboard mcdndeta nsd     $@
	allAboard mcdndeta inelgt0 $@
	allAboard dndeta   inel    $@
	allAboard dndeta   nsd     $@
	allAboard dndeta   inelgt0 $@
    else
	allAboard mcdndeta "" $@
	allAboard dndeta   "" $@
    fi
}
dndetas_upload()
{
    echo "Upload does not make sense for dN/deta jobs"
}
dndetas_draw() 
{
    if test $sys -eq 1 ; then 
	dndeta_draw ${name}_mcdndeta_inel_${now}
	dndeta_draw ${name}_mcdndeta_nsd_${now}
	dndeta_draw ${name}_mcdndeta_inelgt0_${now}
	dndeta_draw ${name}_dndeta_inel_${now} 
	dndeta_draw ${name}_dndeta_nsd_${now} 
	dndeta_draw ${name}_dndeta_inelgt0_${now} 
    else
	dndeta_draw ${name}_mcdndeta_${now}
	dndeta_draw ${name}_dndeta_${now} 
    fi
}

# --- Run all MultDists -------------------------------------------
multdists()
{
    allAboard mcmultdists "" $@
    allAboard multdists   "" $@
}
multdists_upload()
{
    echo "Upload does not make sense for dN/deta jobs"
}
multdists_draw() 
{
    (cd ${name}_mcmultdists_${now} && draw Summarize.C)
    (cd ${name}_multdists_${now}   && draw Summarize.C)
}

# --- Collect PDFs ---------------------------------------------------
collect()
{
    out=${name}_pdfs_${now}
    rm -rf $out
    mkdir -p ${out}
    dirs="corr eloss aod dndeta dndeta_inel dndeta_nsd dndeta_inelgt0 multdists"
    for d in ${dirs} ; do 
	for m in "" "mc" ; do 
	    dir=${name}_${m}${d}_${now}
	    M=
	    case x$m in 
		x)   M=real ;; 
		xmc) M=simu ;;
	    esac
	    if test ! -d $dir ; then 
		# echo "Directory ${dir} doesn't exist"
		continue
	    fi
	    # echo "Will look in $dir"
	    files=
	    case $d in 
		corr)      files="forward_mccorr.pdf" ;; 
		eloss)     files="corrs*.pdf" ;; 
		aod)       files="forward.pdf" ;; 
		dndeta*)   files="forward_dndeta.pdf dNdeta*.pdf" ;; 
		multdists) files="forward_multdists.pdf" ;;
		*) echo "Unknown directory type: $d" > /dev/stder 
		    continue 
		    ;;
	    esac
	    for f in $files ; do 
		ff=$dir/$f
		tgt=
		case $ff in 
		    */forward_mccorr.pdf)    tgt=summary_mccorr.pdf ;; 
		    */forward.pdf)           tgt=summary_${d}_${M}.pdf ;; 
		    */forward_dndeta.pdf)    tgt=summary_${d}_${M}.pdf ;; 
		    */forward_multdists.pdf) tgt=summary_${d}_${M}.pdf ;; 
		    */corr*.pdf)             tgt=summary_${d}_${M}.pdf ;; 
		    */dNdeta*.pdf)           tgt=${d}_${M}.pdf ;;
		    *) echo "Don't know how to deal with $ff" >/dev/stderr 
			continue
			;;
		esac
		# printf "%100s -> %s\n" $ff $tgt
		cp $ff $out/$tgt
	    done
	done
    done 
    (cd ${out} && pdfjoin -q -o tmp.pdf \
	--pdftitle "${name} summary ($now)" \
	--twoside summary_*.pdf dndeta_*.pdf && \
	pdfnup -q --nup 2x1 -o ${name}_summary_${now}.pdf tmp.pdf && \
	rm -f tmp.pdf)
    (cd ${out} && pdfjoin -q -o tmp.pdf \
	--pdftitle "${name} dN/deta ($now)" \
	--twoside dndeta_*.pdf && \
	pdfnup -q --nup 2x1 -o ${name}_dndeta_${now}.pdf tmp.pdf && \
	rm -f tmp.pdf)
    echo "Made ${name}_summary_${now}.pdf and ${name}_dndeta_${now}.pdf"
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
	-r|--run)          run=$opt      ;; 
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
	-W|--workers)      nwrks=${opt}  ;;
	-N|--noact)        noact=1       ;;
	-c|--corrections)  corrs=$opt    ;;
	-a|--par)          par=1         ;;
	-u|--url-opts)     uuopts="$opt" ;;
	-h|--help)         usage         ; exit 0 ;; 
	-H|--manual)       manual        ; exit 0 ;;
	-i|--inel-eff)     inel_eff=$opt ;;
	-0|--inelgt0-eff)  inelgt0_eff=$opt ;;
	-v|--nsd-eff)      nsd_eff=$opt ;;
        --)                break ;;
	*) echo "$0: Unknown option $arg"  ; exit 1 ;; 
    esac
done 

# --- Check settings -------------------------------------------------
check $what

# --- Select what to do ----------------------------------------------
func=
case $what in 
    setup)    setup ; exit 0 ;; 
    clean)    cleanup ; exit 0 ;;
    corr*)    func=corrs;; 
    aod*)     func=aods ;; 
    dndeta*)  func=dndetas ;; 
    multdist*)func=multdists ;;
    collect*) func=collect ;;
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

