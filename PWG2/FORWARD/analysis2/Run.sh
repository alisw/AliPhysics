#!/bin/bash 

ana=$ALICE_ROOT/PWG2/FORWARD/analysis2
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
tit=
dopass1=0
dopass2=0
dopass3=0
pass1=MakeAOD.C
pass2=MakedNdeta.C
pass3=DrawdNdeta.C++g
output1=forward.root
output2=forward_dndeta.root
gdb_script=$ALICE_ROOT/PWG2/FORWARD/analysis2/gdb_cmds
max_rotate=10

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

TYPE is a comma or space separated list of 
 
  INEL	      Inelastic triggers (V0A|V0C|SPD)
  INEL>0      As above + N_ch > 0 in -0.5<eta<+0.5
  NSD         Non-single diffractive ((VOA&VOC)|N_ch > 5 -1.9<eta<+1.9)

If NWORKERS is 0, then the analysis will be run in local mode. 
EOF
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
	-b|--batch)           batch=`toggle $batch`   ;; 
	-P|--proof)           proof=$2	      ; shift ;; 
	-M|--mc)              mc=`toggle $mc`   ;; 
	-g|--gdb)             gdb=`toggle $gdb`   ;; 
	-v|--vz-min)          vzmin=$2         ; shift ;; 
	-V|--vz-max)          vzmax=$2         ; shift ;; 
	-E|--eloss)           pass1=MakeELossFits.C 
	                      pass2=scripts/ExtractELoss.C
	                      pass3=scripts/DrawAnaELoss.C 
	                      output1=forward_eloss.root 
			      dopass2=1 
			     ;;
	-t|--type)           
	    if test "x$type" = "x" ; then type=$2 ; else type="$type|$2"; fi
	    shift ;;
	*) echo "$0: Unknown option '$1'" >> /dev/stderr ; exit 1 ;;
    esac
    shift
done 

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
    rotate AliAOD.root
    rotate ${output1}

    if test $gdb -gt 0 ; then 
	export PROOF_WRAPPERCMD="gdb -batch -x ${gdb_script} --args"
    fi
    echo "Running aliroot ${opts} ${opts1} ${ana}/${pass1}\(\".\",$nev,$proof,$mc\) $redir"
    if test $batch -gt 0 ; then 
	aliroot $opts $opts1 ${ana}/${pass1}\(\".\",$nev,$proof,$mc\) 2>&1 | tee ${base}.log
    else 
	aliroot $opts $opts1 ${ana}/${pass1}\(\".\",$nev,$proof,$mc\)
    fi
    fail=$?
    if  test $fail -gt 0  ; then 
        echo "Return value $fail not 0" ; exit $fail 
    fi
    if test ! -f ${output1} ; then 
	echo "$output1 not made" ; exit 1; 
    fi
    if test ! -f AliAOD.root ; then 
	echo "No AOD creates" ; exit 1;
    fi
    echo "Pass 1 done"
fi

if test $dopass2 -gt 0 ; then
    rotate ${output2}

    args=(\(\".\",$nev,\"$type\",$vzmin,$vzmax,$proof\))
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args=(\(\"${output1}\"\))
    fi
    echo We are Running aliroot ${opts} ${opts1} ${ana}/${pass2}${args}
    aliroot ${opts} ${opts1} ${ana}/${pass2}${args}

    fail=$? 
    if test $fail -gt 0 ; then 
	echo "Return value $fail not 0" ; exit $fail 
    fi
    if test ! -f ${output2} ; then 
	echo "$output2 not made" ; exit 1; 
    fi
    echo "Pass 2 done"
fi

if test $dopass3 -gt 0 ; then
    tit=`echo $tit | tr ' ' '@'` 
    args=(\(\"${output2}\"\))
    if test "x$pass1" = "xMakeELossFits.C" ; then 
	args=(\(\"${output1}\"\))
    fi
    
    echo "Running aliroot ${opts} ${opts1} ${ana}/${pass3}${args}"
    aliroot ${opts} ${ana}/${pass3}${args}
fi				 


#
# EOF
#
