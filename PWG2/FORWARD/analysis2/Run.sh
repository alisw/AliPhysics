#!/bin/bash 

ana=$ALICE_ROOT/PWG2/FORWARD/analysis2
esddir="."
nev=-1
rebin=1
vzmin=-10
vzmax=10
batch=0
gdb=0
proof=0
mc=0
type=INEL
cms=900
hhd=1
comp=1
cent=0
tit=
others=0
published=1
ratios=1
asymm=1 
scheme="full"
dopass1=0
dopass2=0
dopass3=0
pass1=MakeAOD.C
pass2=MakedNdeta.C
pass3=DrawdNdeta.C++g
output1=forward.root
output2=forward_dndeta.root
outputs1="${output1} AliAOD.root event_stat.root EventStat_temp.root"
outputs2="${output2}"
gdb_script=$ALICE_ROOT/PWG2/FORWARD/analysis2/gdb_cmds
max_rotate=10
name=`date +analysis%Y%m%d_%H%M`
pass2dir=./

usage()
{
cat<<EOF
Usage: $0 [OPTIONS]

Do Pass1 and Pass2 on ESD files in current directory.  

Options:
	-h,--help		This help                  
	-n,--events N		Number of events            ($nev)
	-1,--pass1 		Run pass 1, only AOD        ($dopass1)
	-2,--pass2		Run pass 2, only Hists      ($dopass2)
	-3,--pass3		Draw results                ($dopass3)
	-v,--vz-min CM          Minimum value of vz         ($vzmin)
	-V,--vz-max CM          Maximum value of vz         ($vzmax)
	-t,--trigger TYPE       Select trigger TYPE         ($type)
	-b,--batch              Do batch processing         ($batch)
	-P,--proof NWORKERS	Run in PROOF(Lite) mode     ($proof)
	-M,--mc			Run over MC data            ($mc)
	-g,--gdb		Run in GDB mode    	    ($gdb)
	-E,--eloss		Run energy loss script      
        -r,--rebin              Rebin factor                ($rebin)
        -C,--use-centrality     Run centrality task         ($cent)
	-O,--show-older		Show older data	            ($others)
	-J,--show-published	Show ALICE published data   ($published)
	-R,--show-ratios	Show ratios to other data   ($ratios)
	-Z,--show-asymmetry	Show asymmetry 		    ($asymm)
	-S,--scheme SCHEME	Normalisation scheme	    ($scheme)
	-T,--title STRING       Title on plots              ($tit)
	-N,--name STRING        Name of analysis            ($name)
	-I,--input-dir PATH     Path to input ESD data      ($esddir)

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

If NWORKERS is 0, then the analysis will be run in local mode. 
EOF
}

test_ralien()
{
    aliroot -l -b <<EOF > /dev/null 2>&1
int ret = gSystem->Load("libRAliEn");
gApplication->Terminate(ret);
EOF
    ret=$?
    return $ret
}

toggle()
{
    echo $((($1+1)%2))
}

rotate() 
{
    fname=$1 
    if test -f ${fname}.${max_rotate} ; then 
	echo "Maximum number of rotations - $max_rotate - for $fname found" \
	    > /dev/stderr 
	exit 1
    fi
    let max=$max_rotate-1
    for i in `seq $max -1 1` ; do 
	if test ! -f ${fname}.${i} ; then continue ; fi 
	let newn=$i+1
	echo "Moving ${fname}.$i to ${fname}.$newn"
	mv ${fname}.$i ${fname}.$newn
    done
    if test -f $fname ; then 
	echo "Moving ${fname} to ${fname}.1"
	mv $fname ${fname}.1 
    fi
}


while test $# -gt 0 ; do
    case $1 in 
	-h|--help)            usage            ; exit 0;; 
	-n|--events)          nev=$2           ; shift ;; 
	-3|--pass3|-D|--draw) dopass3=`toggle $dopass3`   ;; 
	-2|--pass2|-H|--hist) dopass2=`toggle $dopass2`   ;; 
	-1|--pass1|-A|--aod)  dopass1=`toggle $dopass1`   ;; 
	-b|--batch)           batch=`toggle $batch`       ;; 
	-P|--proof)           proof=$2	          ; shift ;; 
	-C|--use-centrality)  cent=`toggle $cent` ;;
	-M|--mc)              mc=`toggle $mc`     ;; 
	-g|--gdb)             gdb=`toggle $gdb`   ;;
	-r|--rebin)           rebin=$2            ; shift ;;
	-v|--vz-min)          vzmin=$2            ; shift ;; 
	-V|--vz-max)          vzmax=$2            ; shift ;; 
	-E|--eloss)           pass1=MakeELossFits.C 
	                      pass2=scripts/ExtractELoss.C
	                      pass3=scripts/DrawAnaELoss.C 
	                      output1=forward_eloss.root 
			      outputs1="${output1} event_stat.root EventStat_temp.root"
			      outputs2=""
			      dopass2=1 
			      ;;
	-O|--show-older)      others=`toggle $others`	;;
	-J|--show-published)  published=`toggle $published`	;;
	-R|--show-ratios)     ratios=`toggle $ratios`	;;
	-Z|--show-asymmetry)  asymm=`toggle $asymm`	;;
	-S|--scheme)          scheme=`echo $2 | tr ' ' ','` ; shift ;;
	-T|--title)           tit=$2 ; shift ;; 
	-N|--name)            name=$2 ; shift ;; 
	-I|--input-dir)       esddir=$2 ; shift ;;
	-t|--type)           
	    #if test "x$type" = "x" ; then type=$2 ; else type="$type|$2"; fi
	    type=$2
	    shift ;;
	*) echo "$0: Unknown option '$1'" >> /dev/stderr ; exit 1 ;;
    esac
    shift
done 

if test "x$name" != "x" && test_ralien ; then 
    echo "AliEn plug-in available - output will be in $name"
    pass2dir=${name}/
fi

if test $nev -lt 0 ; then 
    base=dndeta_xxxxxxx
else 
    base=`printf dndeta_%07d $nev`
fi
opts="-l -x"
opts1=""
redir=

if test $dopass2 -gt 0 ; then 
    opts1="-q" 
fi
if test $batch -gt 0 ; then 
    opts="-b -q $opts"
    redir="2>&1 | tee ${base}.log"
    echo "redir=$redir"
fi 
if test $dopass1 -gt 0 ; then 
    for i in ${outputs1} ; do 
	rotate ${pass2dir}${i}
    done

    if test $gdb -gt 0 ; then 
	export PROOF_WRAPPERCMD="gdb -batch -x ${gdb_script} --args"
    fi
    echo "Running aliroot $opts $opts1 ${ana}/${pass1}\(\"${esddir}\",$nev,$proof,$mc,$cent,\"${name}\"\)"
    if test $batch -gt 0 ; then 
	aliroot $opts $opts1 ${ana}/${pass1}\(\"${esddir}\",$nev,$proof,$mc,$cent,\"${name}\"\) 2>&1 | tee ${base}.log
    else 
	aliroot $opts $opts1 ${ana}/${pass1}\(\"${esddir}\",$nev,$proof,$mc,$cent,\"${name}\"\)
    fi
    fail=$?
    if  test $fail -gt 0  ; then 
        echo "Return value $fail not 0" ; exit $fail 
    fi
    for i in ${outputs1} ; do 
	if test ! -f ${pass2dir}${i} ; then 
	    echo "File ${i} in ${pass2dir} not generated"
	    exit 1
	fi
	ls -l ${pass2dir}/${i}
    done
    echo "Pass 1 done"
fi

if test $dopass2 -gt 0 ; then
    for i in ${outputs2} ; do 
	rotate ${pass2dir}${i} 
    done 

    args=(\(\"${pass2dir}\",$nev,\"$type\",$cent,\"$scheme\",$vzmin,$vzmax,$proof,\"$name\"\))
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args=(\(\"${pass2dir}${output1}\"\))
    fi
    echo We are Running aliroot ${opts} ${opts1} ${ana}/${pass2}${args}
    aliroot ${opts} ${opts1} ${ana}/${pass2}${args}

    fail=$? 
    if test $fail -gt 0 ; then 
	echo "Return value $fail not 0" ; exit $fail 
    fi
    for i in ${outputs2} ; do 
	if test ! -f ${pass2dir}${i} ; then 
	    echo "File ${i} in ${pass2dir} not generated"
	    exit 1
	fi
	ls -l ${pass2dir}/${i}
    done
    echo "Pass 2 done"
fi

if test $dopass3 -gt 0 ; then
    tit=`echo $tit | tr ' ' '@'` 
    flags=0
    if test $others    -gt 0 ; then let flags=$(($flags|0x1)); fi
    if test $published -gt 0 ; then let flags=$(($flags|0x2)); fi
    if test $ratios    -gt 0 ; then let flags=$(($flags|0x4)); fi
    if test $asymm     -gt 0 ; then let flags=$(($flags|0x8)); fi

    args=(\(\"${pass2dir}${output2}\"\,${flags},\"$tit\",$rebin \))
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args=(\(\"${pass2dir}${output1}\"\))
    fi
    
    echo "Running aliroot ${opts} ${opts1} ${ana}/${pass3}${args}"
    aliroot ${opts} ${ana}/${pass3}${args}
fi				 


#
# EOF
#
