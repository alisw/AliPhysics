#!/bin/bash 

# General options 
batch=0
proof=0
cent=0
dopass1=0
dopass2=0
dopass3=0
max_rotate=10
prog=aliroot

# Input (data and scripts)
ana=$ALICE_ROOT/PWG2/FORWARD/analysis2
esddir="."
nev=-1

# Pass 1 options 
mc=0
scheme="full"

# Pass 2 event selection
vzmin=-10
vzmax=10
type=INEL
mcfilename="none"

# Pass 3 visualisation options 
rebin=1
tit=
others=
ratios=1
asymm=1 
rings=0
syserr=1

# Derived variables 
pass1=MakeAOD.C
pass2=MakedNdeta.C
pass3=DrawdNdeta.C++
output1=forward.root
output2=forward_dndeta.root
outputs1="${output1} AliAOD.root event_stat.root EventStat_temp.root"
outputs2="${output2}"
gdb_script=$ALICE_ROOT/PWG2/FORWARD/analysis2/gdb_cmds
name=`date +analysis%Y%m%d_%H%M`
pass2dir=./

#_____________________________________________________________________
# Print usage
usage()
{
cat<<EOF
Usage: $0 [OPTIONS]

Do Pass1 and Pass2 on ESD files in current directory.  

Options:
   General options:
	-h,--help		This help                  
	-b,--batch              Do batch processing         ($batch)
	-P,--proof NWORKERS	Run in PROOF(Lite) mode     ($proof)
	-1,--pass1 		Run pass 1 (AOD)            ($dopass1)
	-2,--pass2		Run pass 2 (Hists)          ($dopass2)
	-3,--pass3		Run pass 3 (Vizualisation)  ($dopass3)
	-n,--events N		Number of events            ($nev)
	-g,--gdb		Run in GDB mode    	    ($gdb)
	-N,--name STRING        Name of analysis            ($name)
	-E,--eloss		Run energy loss script      
	-A,--script-dir		Script directory            ($ana)
	-a,--program		Program to use		    ($prog)
   Pass 1 options and event selection
        -C,--use-centrality     Run centrality task         ($cent)
	-M,--mc			Run over MC data            ($mc)
	-I,--input-dir PATH     Path to input ESD data      ($esddir)
   Pass 2 options 
	-v,--vz-min CM          Minimum value of vz         ($vzmin)
	-V,--vz-max CM          Maximum value of vz         ($vzmax)
	-t,--trigger TYPE       Select trigger TYPE         ($type)
	-S,--scheme SCHEME	Normalisation scheme	    ($scheme)
	-q,--mcfilename		Result of MC analysis	    ($mcfilename)
   Pass 3 options 
        -r,--rebin              Rebin factor                ($rebin)
	-O,--others WHICH	Show other data	            ($others)
	-J,--show-rings		Show individual rings       ($rings)
	-R,--show-ratios	Show ratios to other data   ($ratios)
	-Z,--show-asymmetry	Show asymmetry 		    ($asymm)
	-G,--show-syserror      Show systematic errors      ($syserr)
	-T,--title STRING       Title on plots              ($tit)

TYPE is a comma or space separated list of 
 
  INEL	      Inelastic triggers (V0A|V0C|SPD)
  INEL>0      As above + N_ch > 0 in -0.5<eta<+0.5
  NSD         Non-single diffractive ((VOA&VOC)|N_ch > 5 -1.9<eta<+1.9)

SCHEME is a comma or space separated list of 

  NONE          No event-level normalization except trivial one 
  EVENT         Event-level normalization 
  BACKGROUND    Not implemented yet 
  SHAPE         Shape correction 
  TRIGGER       Trigger efficiency 
  FULL          Same as EVENTLEVEL,BACKGROUND,SHAPE,TRIGGER

WHICH is a comma or space separated list of 

  UA5		UA5 data (900GeV, ppbar, INEL, NSD)
  CMS           CMS data 
  ALICE         ALICE published data
  PYTHIA        Event generator data 

If NWORKERS is 0, then the analysis will be run in local mode. 
EOF
}

#_____________________________________________________________________
test_ralien()
{
    aliroot -l -b <<EOF > /dev/null 2>&1
int ret = gSystem->Load("libRAliEn");
gApplication->Terminate(ret);
EOF
    ret=$?
    return $ret
}

#_____________________________________________________________________
# Toggle a value 
toggle()
{
    echo $((($1+1)%2))
}

#_____________________________________________________________________
# Rotate files 
rotate() 
{
    fname=$1 
    if test -f ${fname}.${max_rotate} ; then 
	# echo "Removing ${fname}.${max_rotate}"
	rm -f ${fname}.${max_rotate}
	# echo "Maximum number of rotations - $max_rotate - for $fname found" \
	#     > /dev/stderr 
	# exit 1
    fi
    let max=$max_rotate-1
    for i in `seq $max -1 1` ; do 
	if test ! -f ${fname}.${i} ; then continue ; fi 
	let newn=$i+1
	# echo "Moving ${fname}.$i to ${fname}.$newn"
	mv ${fname}.$i ${fname}.$newn
    done
    if test -f $fname ; then 
	# echo "Moving ${fname} to ${fname}.1"
	mv $fname ${fname}.1 
    fi
}

#_____________________________________________________________________
#
# Function to run pass 
#
# Arguments (in order)
#    isBatch:     Should we do batch processing 
#    output:      Main output file 
#    outputs:     All outputs 
#    outdir:      Where the output will be put 
#    notlast:     Not last pass 
#    script:      Script name 
#    args:        Arguments for script 
#
run_pass()
{
    isbatch=$1    ; shift 
    output=$1     ; shift
    outputs=$1    ; shift 
    outdir=$1     ; shift 
    notLast=$1    ; shift 
    script=$1     ; shift 
    args=$1       

    # --- Log file name ----------------------------------------------
    if test "x$output" = "x" ; then 
	log=analysis.log
    else
	log=`dirname $output`/`basename $output .root`.log
    fi

    # --- Make options for AliROOT -----------------------------------
    opts=
    if test $isbatch -gt 0 || test $notLast -gt 0 ; then 
	opts="-q"
    fi
    if test $isbatch -gt 0 ; then 
	opts="-b $opts"
    fi 

    # --- Rotate output file -----------------------------------------
    for i in ${outputs} ${log} ; do 
	rotate ${outdir}${i}
    done

    # --- Some print out ---------------------------------------------
    cat <<-EOF
	Pass parameters: 
	 Batch mode:          	$isBatch 
	 Main output:	     	$output
	 All outputs:		$outputs
	 Output directory:	$outdir
	 More to do:		$notLast
	 Script to run:		$script 
	 Script arguments:	$args
	 Log file:              $log
	 AliROOT options:	$opts
	EOF

    # --- Run AliROOT ------------------------------------------------
    echo "Running ${prog} $opts ${script}${args}"
    if test $isbatch -gt 0 || test $notLast -gt 0 ; then 
	${prog} ${opts} ${script}"${args}" 2>&1 | tee ${log}
    else 
	${prog} ${opts} ${script}"${args}"
    fi
    fail=$?

    # --- Check exit conditions --------------------------------------
    if  test $fail -gt 0  ; then 
        echo "Returned $fail" 
	exit $fail 
    fi
    for i in ${outputs} ; do 
	if test ! -f ${outdir}${i} ; then 
	    echo "File ${i} not generated in ${outdir}"
	    exit 1
	fi
    done
    echo "Success (log in $log)"
}

#_____________________________________________________________________
# Loop over arguments 
while test $# -gt 0 ; do
    case $1 in 
	# General options 
	-h|--help)            usage  ; exit 0		  ;; 
	-b|--batch)           batch=`toggle $batch`       ;; 
	-P|--proof)           proof=$2	; shift 	  ;; 
	-1|--pass1|--aod)     dopass1=`toggle $dopass1`   ;; 
	-2|--pass2|--hist)    dopass2=`toggle $dopass2`   ;; 
	-3|--pass3|--draw)    dopass3=`toggle $dopass3`   ;; 
	-n|--events)          nev=$2 ; shift 		  ;;
	-g|--gdb)             gdb=`toggle $gdb`   	  ;;
	-N|--name)            name=$2 ; shift 		  ;; 
	-A|--script-dir)      ana=$2; shift 		  ;;
	-a|--program)	      prog=$2 ; shift 		  ;;
	-E|--eloss)           pass1=MakeELossFits.C 
	                      pass2=scripts/ExtractELoss.C
	                      pass3=scripts/DrawAnaELoss.C 
	                      output1=forward_eloss.root 
			      outputs1="${output1} event_stat.root EventStat_temp.root"
			      outputs2=""
			      dopass2=1 
			      ;;
	# Pass 1 options
	-C|--use-centrality)  cent=`toggle $cent` 	;;
	-M|--mc)              mc=`toggle $mc`     	;; 
	-I|--input-dir)       esddir=$2 ; shift 	;;
	# Pass 2 options 
	-v|--vz-min)          vzmin=$2 ; shift 		;; 
	-V|--vz-max)          vzmax=$2 ; shift 		;; 
	-S|--scheme)          scheme=`echo $2 | tr ' ' ','` ; shift ;;
	-t|--type)            type=$2  ; shift ;;
	-q|--mcname)          mcfilename=$2 ; shift ;;
	# PAss 3 options 
	-r|--rebin)           rebin=$2  ; shift 	;;
	-O|--others)          others=`echo $2 | tr ',' ' '` ; shift	;;
	-R|--show-ratios)     ratios=`toggle $ratios`	;;
	-J|--show-rings)      rings=`toggle $rings`	;;
	-Z|--show-asymmetry)  asymm=`toggle $asymm`	;;
	-G|--show-syserror)   syserr=`toggle $syserr`	;;
	-T|--title)           tit=$2 ; shift ;;
	*) echo "$0: Unknown option '$1'" >> /dev/stderr ; exit 1 ;;
    esac
    shift
done 

#_____________________________________________________________________
# Check for RAlien if needed
if test "x$name" != "x" && test_ralien ; then 
    echo "AliEn plug-in available - output will be in $name"
    pass2dir=${name}/
fi

#_____________________________________________________________________
# Pass 1 
if test $dopass1 -gt 0 ; then 
    args="(\"${esddir}\",$nev,$proof,$mc,$cent,\"${name}\")"
    echo "Args=$args"
    run_pass ${batch} ${output1} "${outputs1}" "${pass2dir}" ${dopass2} \
	${ana}/${pass1} ${args}
    echo "Pass 1 done"
fi

#_____________________________________________________________________
# Pass 2 
if test $dopass2 -gt 0 ; then
    args="(\"${pass2dir}\",$nev,\"$type\",$cent,\"$scheme\",$vzmin,$vzmax,$proof,\"$name\",\"$mcfilename\")"
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args=(\(\"${pass2dir}${output1}\"\))
    fi
    
    run_pass ${batch} ${output2} "${outputs2}" "${pass2dir}" ${dopass3} \
	${ana}/${pass2} ${args}
    echo "Pass 2 done"
fi

#_____________________________________________________________________
# Pass 3 
if test $dopass3 -gt 0 ; then
    tit=`echo $tit | tr ' ' '@'` 
    flags=0
    if test $ratios    -gt 0 ; then let flags=$(($flags|0x1)); fi
    if test $asymm     -gt 0 ; then let flags=$(($flags|0x2)); fi
    if test $syserr    -gt 0 ; then let flags=$(($flags|0x4)); fi
    if test $rings     -gt 0 ; then let flags=$(($flags|0x8)); fi
    which=0
    others=`echo $others | tr ',' ' ' | tr '[a-z]' '[A-Z]'` 
    for i in $others ; do 
	case $i in 
	    UA5)	let which=$(($which|0x1)) ;; 
	    CMS)	let which=$(($which|0x2)) ;; 
	    ALICE)	let which=$(($which|0x4)) ;; 
	    PYTHIA|MC)	let which=$(($which|0x8)) ;; 
	esac
    done 

    args="(\"${pass2dir}${output2}\",\"$tit\",$rebin,${which},${flags})"
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args="(\"${pass2dir}${output1}\")"
    fi

    run_pass ${batch} "" "" "" 0 ${ana}/${pass3} ${args}
fi				 


#
# EOF
#
