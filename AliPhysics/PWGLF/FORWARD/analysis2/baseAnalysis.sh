#
# This file contain common functions used by liteAnalysis.sh,
# gridAnalysis.sh, etc. 
#

# === Variables ======================================================
fwd_dir=$ALICE_PHYSICS/PWGLF/FORWARD/analysis2
dotconf=.config
here=${PWD}
name=
now=
sys=0
snn=0
field=0
corrs=
noact=0
uuopts=
# Latest: 0.91
inel_eff=1
# Latest: 0.94
nsd_eff=1
inelgt0_eff=1

# === Misc. ==========================================================
dummy()
{
    local f=$1 
    echo "$f: Dummy function. Please implement in $0" > /dev/stderr 
    exit 1
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

# === Utilities to execute scripts ===================================
# --- Run script -----------------------------------------------------
script()
{
    local scr=$1 ; shift 
    local args=$1 ; shift
    echo "Will run aliroot -l -b -q $scr($args)"
    if test $noact -gt 0 ; then return ; fi
    aliroot -l -b <<-EOF
	.x $scr($args)
	.q
	EOF
}
# --- Extract corrections --------------------------------------------
terminate()
{
    echo "Nothing to do for terminate"  
}
 # --- Post processing -----------------------------------------------
post()
{
    ./post.sh $@ 
}
# --- Extract corrections --------------------------------------------
_extract()
{
    local scr=${1:-Extract.C}
    if test ! -f ${scr} ; then 
	scr="../${scr}" 
	if test ! -f ${scr} ; then 
	    echo "Extract script not found in `pwd` or parent" > /dev/stderr 
	    exit 1
	fi
    fi
    if test -f .extract ; then 
	echo "Aldready extracted in `basename $PWD`" 
	return 0 
    fi
    echo "= Extracting"
    script $scr > /dev/null 2>&1
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
    _extract
    upload
}
# --- Draw -----------------------------------------------------------
_draw()
{
    script $@
}
# --- Generic draw ---------------------------------------------------
draw()
{
    _draw $@ 
}
# --- Draw dN/deta results -------------------------------------------
_dndeta_draw()
{
    local d=$1 ; shift
    local scr=${1:-Draw.C}
    echo "=== $d ================================================"
    (cd $d && \
	draw ${fwd_dir}/DrawdNdetaSummary.C && \
	draw ${scr})
    echo "Back to `pwd`"
}
# --- Draw dN/deta ---------------------------------------------------
dndeta_draw()
{
    _dndeta_draw $@ 
}
# === Task functions =================================================
# --- Common options help --------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 --what=OPERATION [--step=STEP] [OPTIONS]

General options:
  -s,--step=STEP            Run stage ($step)
  -w,--what=TRAINS          What to do 
  -M,--man                  Show the manual  
  -N,--noact                Show what will be done 

Options for 'setup' operation:
  -n,--name=STRING          Base name of jobs ($name)
  -S,--sys=SYSTEM           Collision system ($sys)
  -E,--snn=ENERGY           Center of mass energy per nuclean pair ($snn)
  -F,--field=FIELD          L3 magnetic field ($field)
  -c,--corrections=DIR      Directory where corrections are stored ($corrs)
  -u,--url-opts=OPTIONS	    Additional user options ($uuopts)
  -i,--inel-eff=EFF         Set INEL efficiency (only pp - $inel_eff)
  -0,--inelgt0-eff=EFF      Set INEL>0 efficiency (only pp - $inelgt0_eff)
  -v,--nsd-eff=EFF          Set NSD (to V0-AND) efficiency (only pp $nsd_eff)
EOF
    setup_usage
    cat <<EOF

TRAINS is one of

  clean       Clean directory
  setup       Do intial setup 
  corrs       Generate corrections 
  aods        Generate AODs 
  dndeta      Generate dNdeta 
  multdists   Generate P(Nch)
  flow        Generate v_n{m} - not implemented yet 

and must be executed in that order.  STEP is one of 

  full        Run the analysis 
  terminate   Terminate the job (may need iterations)
  upload      Upload corrections (only for TRAINS=corrs)
  draw        Draw (partial) results (not for TRAINS=corrs)
EOF
    
}
setup_usage()
{
    dummy setup_usage
}

# === Setup functions ================================================
# --- Common setup code ----------------------------------------------
_setup()
{
    local lhandled=yes
    declare -a larg=$@
    while test $# -gt 0 ; do 
	arg=$1 
	opt=
	case $1 in 
	    --*=*) 
		arg=`echo $1 | sed 's/=.*//'` ;
		opt=`echo $1 | sed 's/--[^=][^=]*=//'` 
		;;
	    --*)
		;;
	    -*) opt=$2 ; shift ;; 
	esac
	shift 

	lhandled=yes
	# echo "Base: Processing '$arg' ('$opt')"
	case $arg in 
	    -n|--name)         name=$opt     ;; 
	    -c|--corrections)  corrs=$opt    ;;
	    -i|--inel-eff)     inel_eff=$opt ;;
	    -0|--inelgt0-eff)  inelgt0_eff=$opt ;;
	    -v|--nsd-eff)      nsd_eff=$opt ;;
	    -u|--url-opts)     uuopts="$opt" ;;
	    -S|--sys)          sys=`echo $opt | tr '[A-Z]' '[a-z]'` ;;
	    -E|--snn)          snn=$opt      ;; 
	    -F|--field)        field=$opt    ;;
	    *) lhandled=no ;;
	esac
	if test x$lhandled = "xno" ; then 
	    handle_setup_option "$arg" "$opt" 
	fi
    done

    check setup

    # Set the date/time string 
    now=`date '+%Y%m%d_%H%M'` 

    # define our outputs
    echo "=== Define outputs"
    outputs

    # Dump configuration to file 
    echo "=== Dump configuration"
    dump_conf $larg

    # Set-up for corrections 
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
	  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
	  if (!gROOT->GetClass("AliOADBForward"))
	    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
	  gROOT->LoadMacro(Form("%s/corrs/ForwardOADBGui.C", fwd));
	  
	  AliOADBForward* db = new AliOADBForward;
	  db->Open("fmd_corrections.root", "*");
	  
	  ForwardOADBGui(db);

	  return db;
	}
	EOF
	echo "=== Make acceptance corrections" 
	(cd ${name}_acc_${now} && \
	    accGen `run_for_acc` && \
	    upload )
   fi
   for i in fmd_corrections.root spd_corrections.root deadstrips.C ; do 
       if test ! -f ${corrdir}/$i ; then continue ; fi 
       echo "Linking ${corrdir}/$i here"
       ln -fs ${corrdir}/$i . 
   done
   print 

    
}
# --- Default implementation -----------------------------------------
setup()
{
    echo "Default implementation"
    _setup $@ 
}
# --- dummy handler of setup options ---------------------------------
handle_setup_option()
{
    dummy handle_setup_option
}
# --- Dump configuration to file -------------------------------------
dump_conf()
{
   cat > ${dotconf} <<-EOF
	# Generated by command
	# 
	EOF
   echo -n "#  $0 " >> $dotconf
   while test $# -gt 0 ; do 
       echo -en " \\" >> $dotconf
       case $1 in 
	   --*) echo -en "\n#\t$1"    >> $dotconf;; 
	   -*)  echo -en "\n#\t$1 $2" >> $dotconf; shift ;; 
	   *)   echo -en "\n#\t$1"    >> $dotconf;;
       esac
       shift 
   done
   cat >> ${dotconf} <<-EOF
	
	# Settings:
	name="$name"
	now=${now}
	# Collision system and similar
	sys=$sys
	snn=$snn
	field=$field
	# Additional URI options 
	uuopts="${uuopts}"
	# Trigger efficiencies - edit here to set them 
	inel_eff=$inel_eff
	inelgt0_eff=$inelgt0_eff
	nsd_eff=$nsd_eff
	EOF
   dump_setup ${dotconf}
   echo "# EOF" >> ${dotconf}
}
# --- Get the run number to use for acceptance -----------------------
run_for_acc()
{
    dummy run_for_acc
}
# --- Dump settings to output file -----------------------------------
dump_setup()
{
    dummy dump_setup
}

# --- Run acceptance generation --------------------------------------
accGen()
{
    check_token
    local run=$1 
    script ${fwd_dir}/corrs/ExtractAcceptance.C "${run}"
    if test -f deadstrips.C ; then 
	cp ${here}/${name}_corrs_${now}/
    fi
}

# === Check function =================================================
check()
{
    local w=$1

    if test "X$name" = X ; then 
	echo "No name specified" > /dev/stderr 
	exit 1
    fi
    if test "X$w" != "Xsetup" && test "x$now" = "x" ; then 
	echo "No date/time specified" > /dev/stderr 
	exit 1
    fi

    check_setup $w

    # sys==0 is OK - autoselect
    case x$sys in 
        xpp|xp-p)              sys=1 ;; 
	xpbpb|xpb-pb|xaa|xa-a) sys=2 ;; 
	xppb|xp-pb|xpa|xp-a)   sys=3 ;;
	x0|x1|x2|x3)                 ;; 
	x)                     sys=0 ;;
	*) echo "$0: Unknown system: $sys" ; exit 1 ;;
    esac

}
check_setup()
{
    dummy check_setup
}

# --- Show the setup -------------------------------------------------
print()
{
    cat <<-EOF
	Name:			${name}
	Collision system:	${sys}
	sqrt(s_NN):		${snn}GeV
	L3 Field:		${field}kG
	Date & time:            ${now}
	Additional URL options: ${uuopts}
	Trigger efficiencies:   
	  INEL:                 ${inel_eff}
	  INEL>0:               ${inelgt0_eff}
	  NSD:                  ${nsd_eff}
	EOF
    print_setup
}
print_setup()
{
    dummy print_setup
}

# === Task functions =================================================
# Modifies 'opts' 
train_opts()
{
    local mc=$1   ; shift
    local type=$1 ; shift 
    local trig=$1 ; shift 

    nme=${name}_${type}

    case $type in 
	*corr)   cl=MakeMCCorrTrain ; mc=1 ;;
	*eloss)  cl=MakeFMDELossTrain ;;  
	*aod)    cl=MakeAODTrain    ; opts="${opts} --corr=." ;;
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
	    if test $mc -gt 0 ; then 
		opts="$opts --mc"
	    fi
	    ;;
	*multdists) 
	    cl=MakeMultDistsTrain 
	    tree=aodTree 
	    pat="*/AliAOD.root"
	    ;;
	*) echo "$0: Unknown type of train: $type" > /dev/stderr ; exit 1 ;;
    esac
    # add centrality flag if we do not know what collision system we're 
    # looking at, or it's PbPb or pPb. 
    case $sys in 
	0|2|3) opts="$opts --cent" ;; 
	1)                         ;;
    esac
    opts="${opts} --include=$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/trains"
    opts="${opts} --date=${now} --class=$cl --name=$nme --verbose=0"
}

# --- Run the train --------------------------------------------------
_allAboard()
{
    local type=$1 ; shift 
    local trig=$1 ; shift 

    local mc=0
    case $type in 
	mc*) mc=1  ;; 
    esac

    opts=
    url=
    train_opts $mc $type $trig 
    url_opts $mc $type $trig

    
    echo "=== Running train: runTrain ${opts} --url=${url} $@" 
    if test $noact -gt 0 ; then return ; fi

    runTrain ${opts} --overwrite --url=${url} $@ 
}
url_opts()
{
    dummy url_opts
}
allAboard()
{
    local what=$1 ; shift 
    local trig=$1 ; shift 
    echo "=== Running _allAboard '$what' '$trig' $@" 
    _allAboard "$what" "$trig" $@
}

# --- Collect PDFs ---------------------------------------------------
collect()
{
    local out=${name}_pdfs_${now}
    rm -rf $out
    mkdir -p ${out}
    dirs="corr eloss aod dndeta dndeta_inel dndeta_nsd dndeta_inelgt0 multdists"
    # Now loop on dirs
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
	    local files=
	    case $d in 
		corr)      files="forward_mccorr.pdf" ;; 
		eloss)     files="forward_elossfits.pdf" ;; 
		aod)       files="forward.pdf" ;; 
		dndeta*)   files="forward_dndeta.pdf dNdeta*.pdf" ;; 
		multdists) files="forward_multdists.pdf" ;;
		*) echo "Unknown directory type: $d" > /dev/stder 
		    continue 
		    ;;
	    esac
	    collect_files "$dir" "$d" "$M" "$out" "$run" "$files"

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
    rm -f last_${name}_pdfs
    ln -s ${out} last_${name}_pdfs
    (cd ${out} && _pres_collect)
}

_pres_collect()
{
    local out=results.tex
    rm -rf $out
    local A=`getent passwd $USER | cut -f5 -d: `
    local D=`echo $now | sed 's,\(....\)\(..\)\(..\)_\(..\)\(..\),\1/\2/\3 \4:\5,'`
    cat <<-EOF > $out
	\\documentclass[compress]{beamer}
	\\usetheme{alice}
	\\usepackage[english,british]{babel}
	\\mode<presentation>
	\\title{Job ${name}}
	\\author{$A} 
	\\date{$D}
	\\begin{document}
	\\aliceTitlePage{}
	EOF
    dirs="dndeta dndeta_inel dndeta_nsd dndeta_inelgt0 multdists"
    for i in ${dirs} ; do 
	l=`ls ${i}_real_*.pdf 2>/dev/null` 
	for r in $l ; do 
	    m=`echo $r | sed 's/real/simu/'` 
	    R=`basename $r .pdf | sed "s/${i}_real_0*//"` 
	    t=`basename $r .pdf | sed 's/dndeta_\(.*\)_real.*/\1/'` 
	    case x$t in 
		x) T="All";;
		xinel)    T="INEL" ;; 
		xinelgt0) T="INEL\\textgreater0" ;; 
		xnsd)     T="NSD" ;; 
		x*)       T="Unknown";;
	    esac
	    if test -f $r && test -f $m ; then 
		cat <<-EOF >> $out
		\\begin{frame}{Run $R --- $T} 
		  \\begin{columns}
		    \\begin{column}{.48\\linewidth}
		      \\textbf{Real}\\\\
		      \\includegraphics[keepaspectratio,width=\\linewidth]{%
		   	$r}
		    \\end{column}
		    \\begin{column}{.48\\linewidth}
		      \\textbf{Simulated}\\\\
		      \\includegraphics[keepaspectratio,width=\\linewidth]{%
		   	$m}
		    \\end{column}
		  \\end{columns}
		\\end{frame}
		EOF
	    fi
	done
    done
    cat <<-EOF >> $out
	\\end{document}
	EOF
    pdflatex $out
}

_collect_files()
{
    local dir=$1 ; shift
    local d=$1 ; shift 
    local M=$1 ; shift 
    local out=$1 ; shift
    local r=$1 ; shift
    # local files="$1"
    for f in $@ ; do 
	ff=$dir/$f
	tgt=`collect_name "$ff" "$d" "$M"`
	if test "x$tgt" = "x" ; then 
	    continue 
	fi
	tgt=`printf "%s_%09d.pdf" $tgt $run` 
	# printf "%100s -> %s\n" $ff $tgt
	if test ! -f $ff ; then 
	    echo "$ff not found - ignored"
	    continue
	fi
	cp $ff $out/$tgt
    done # for f in files
}

collect_name()
{
    local tgt=
    local ff="$1" ; shift
    local d="$1" ; shift 
    local M=$1 ; shift
    case $ff in 
	*/forward_mccorr.pdf)    tgt=summary_mccorr ;; 
	*/forward.pdf)           tgt=summary_${d}_${M} ;; 
	*/forward_dndeta.pdf)    tgt=summary_${d}_${M} ;; 
	*/forward_multdists.pdf) tgt=summary_${d}_${M} ;; 
	*/forward_elossfits.pdf) tgt=summary_${d}_${M} ;; 
	*/dNdeta*.pdf)           tgt=${d}_${M} ;;
	*) echo "Don't know how to deal with $ff" >/dev/stderr 
    esac
    echo $tgt
}

collect_files()
{
    _collect_files $@
}

# --- Clean up -------------------------------------------------------
_cleanup()
{
    for i in acc aod corrs dndeta dndeta_inel dndeta_inelgt0 dndeta_nsd \
	eloss multdist ; do 
	for j in "" mc ; do 
	    if test ! -d ${name}_${j}${i}_${now} ; then 
		continue 
	    fi
	    rm -rf ${name}_${j}${i}_${now} 
	    rm -f last_${name}_${j}${i} 
	done
    done
    rm -f build.log 
    for i in fmd_corrections.root spd_corrections.root deadstrips.C ; do
	if test ! -L $i ; then 
	    continue
	fi
	rm -f fmd_corrections.root
    done
}
cleanup()
{
    _cleanup
}
# === Wrappers =======================================================
# --- Run all correction jobs ----------------------------------------
# This assumes the function allAboard 
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
# Run post processing 
corrs_post()
{
    (cd ${name}_mccorr_${now}  && post ../${name}_corrs_${now})
    (cd ${name}_mceloss_${now} && post ../${name}_corrs_${now} )
    (cd ${name}_eloss_${now}   && post ../${name}_corrs_${now} )
    rm -f fmd_corrections.root spd_corrections.root 
    ln -s ${name}_corrs_${now}/fmd_corrections.root .
    ln -s ${name}_corrs_${now}/spd_corrections.root .
}

# This assumes the function extract_upload 
corrs_upload() 
{
    corrs_post
}
corrs_draw()
{
    corrs_post
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
    aods_post
}
aods_draw() 
{
    aods_post
}
aods_post()
{
    (cd ${name}_mcaod_${now} && post)
    (cd ${name}_aod_${now}   && post)
}
# --- Run all dN/deta jobs -------------------------------------------
dndetas()
{
    if test $sys -eq 1 ; then 
	allAboard mcdndeta inel $@
	allAboard mcdndeta nsd $@
	allAboard mcdndeta inelgt0 $@
	allAboard dndeta inel $@
	allAboard dndeta nsd $@
	allAboard dndeta inelgt0 $@
    else
	allAboard mcdndeta "" $@
	allAboard dndeta   "" $@
    fi
}
dndetas_terminate() 
{
    if test $sys -eq 1 ; then 
	(cd ${name}_mcdndeta_inel_${now}    && terminate)
	(cd ${name}_mcdndeta_nsd_${now}     && terminate)
	(cd ${name}_mcdndeta_inelgt0_${now} && terminate)
	(cd ${name}_dndeta_inel_${now}      && terminate)
	(cd ${name}_dndeta_nsd_${now}       && terminate)
	(cd ${name}_dndeta_inelgt0_${now}   && terminate)
    else
	(cd ${name}_mcdndeta_${now} && terminate)
	(cd ${name}_dndeta_${now}   && terminate)
    fi
}
dndetas_upload()
{
    dndeta_post
}
dndetas_draw() 
{
    dndeta_post
}
dndetas_post() 
{
    if test $sys -eq 1 ; then 
	(cd ${name}_mcdndeta_inel_${now}    && post)
	(cd ${name}_mcdndeta_nsd_${now}     && post)
	(cd ${name}_mcdndeta_inelgt0_${now} && post)
	(cd ${name}_dndeta_inel_${now}      && post)
	(cd ${name}_dndeta_nsd_${now}       && post)
	(cd ${name}_dndeta_inelgt0_${now}   && post)
    else
	(cd ${name}_mcdndeta_${now} && post)
	(cd ${name}_dndeta_${now}   && post)
    fi
}
# --- Run all MultDists -------------------------------------------
multdists()
{
    allAboard mcmultdists "" $@
    allAboard multdists   "" $@
}
multdists_terminate() 
{
    (cd ${name}_mcmultdists_${now} && terminate)
    (cd ${name}_multdists_${now}   && terminate)
}
multdists_upload()
{
    multdists_posts
}
multdists_post() 
{
    (cd ${name}_mcmultdists_${now} && post)
    (cd ${name}_multdists_${now}   && post ../${name}_mcmultdists_${now})
}

# === Driver code ====================================================
# assumes functions:
#  
#  handle_option ARG [OPT]
#  check [WHAT]
#  print_setup
#  setup
#  cleanup
#
runIt()
{
    # --- Source settings if found -----------------------------------
    if test -f $dotconf ; then 
	source $dotconf 
    fi


    # --- Process command line -------------------------------------------
    what=
    step=
    declare -a parg
    local iarg=0
    while test $# -gt 0 ; do
	local arg=$1 
	local opt=
	local handled_arg=1
	local handled_opt=0
	case $1 in 
	    --) shift ; break ;;
	    --*=*) 
		arg=`echo $1 | sed 's/=.*//'` ;
		opt=`echo $1 | sed 's/--[^=][^=]*=//'` 
		;;
	    --*)
		;;
	    -h|-N|-H|-a) ;;
	    -*) opt=$2 ; handled_opt=1 ;; 
	esac

	# echo "Run:  Processing '$arg' ('$opt')"
	case $arg in 
	    -w|--what)         what=`echo $opt | tr '[A-Z]' '[a-z]'`  ;; 
	    -s|--step)         step=`echo $opt | tr '[A-Z]' '[a-z]'`  ;; 
	    -N|--noact)        noact=1       ;;
	    -h|--help)         usage         ; exit 0 ;; 
	    -H|--manual)       manual        ; exit 0 ;;
	    *) handled_arg=0 ;;
	esac
	if test $handled_arg -lt 1 ; then 
	    parg[$iarg]="$1"
	    let iarg=$iarg+1
	    # echo "Pushed $1 -> ${parg[@]}"
	    if test $handled_opt -gt 0 ; then 
		parg[$iarg]="$2"
		let iarg=$iarg+1
		# echo "Pushed $2, gobble $2 -> ${parg[@]}"
		shift
	    fi
	else
	    if test $handled_opt -gt 0 ; then 
		shift 
		# echo "Gobble $1"
	    fi
	fi
	shift 
    done 
    # For backward compatibility 
    while test $# -gt 0 ; do
	parg[$iarg]=$1 
	shift 
	let iarg=$iarg+1
    done
    # --- Select what to do ------------------------------------------
    func=
    case $what in 
	setup)    setup ${parg[@]} ; exit 0 ;; 
	clean)    cleanup ; exit 0 ;;
	corr*)    func=corrs;; 
	aod*)     func=aods ;; 
	dndeta*)  func=dndetas ;; 
	multdist*)func=multdists ;;
	collect*) func=collect ;;    
	*) echo "$0: Unknown operation: $what" > /dev/stderr ; exit 1 ;;
    esac

    # --- Check settings ---------------------------------------------
    check $what

    print

    case $what in 
	setup|clean|collect) step=full ;; 
    esac

    case x$step in 
	x|xfull) ;; 
	xterm*) func=${func}_terminate ;; 
	xup*)   func=${func}_upload ;; 
	xdr*)   func=${func}_draw ;;
	xpo*)   func=${func}_post ;;
	*) echo "$0: Unknown step $step" > /dev/stderr ; exit 1 ;;
    esac
    
    echo "Will execute $func (${parg[@]})" 
    $func ${parg[@]}
}

# 
# EOF
#

