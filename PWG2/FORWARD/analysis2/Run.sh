#!/bin/bash 

ana=$ALICE_ROOT/PWG2/FORWARD/analysis2
nev=-1
noanal=0
nodraw=0
rebin=1
vzmin=-10
vzmax=10
batch=0
gdb=0
proof=0
mc=0
type=
cms=900
hhd=1
comp=1
tit=
pass1=Pass1.C
output=forward.root

usage()
{
cat<<EOF
Usage: $0 [OPTIONS]

Do Pass1 and Pass2 on ESD files in current directory.  

Options:
	-h,--help		This help                  
	-n,--events N		Number of events           ($nev)
	-1,--pass1 		Run only pass 1, no draw   ($nodraw)
	-2,--pass2		Run only pass 2, just draw ($noanal)
	-r,--rebin N		Rebin by N                 ($rebin)
	-v,--vz-min CM          Minimum value of vz        ($vzmin)
	-V,--vz-max CM          Maximum value of vz        ($vzmax)
	-t,--trigger TYPE       Select trigger TYPE        ($type)
	-e,--energy CMS         Center of mass energy      ($cms)
	-b,--batch              Do batch processing        ($batch)
	-P,--proof NWORKERS	Run in PROOF(Lite) mode    ($proof)
	-M,--mc			Run over MC data           ($mc)
	-S,--title STRING       Set the title string       ($tit)
	-g,--gdb		Run in GDB mode    	   ($gdb)
	-H,--hhd		Do comparison to HHD	   ($hhd)
	-O,--other		Do comparison to other	   ($comp)
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


while test $# -gt 0 ; do
    case $1 in 
	-h|--help)           usage            ; exit 0;; 
	-n|--events)         nev=$2           ; shift ;; 
	-2|--pass2)          noanal=`toggle $noanal`   ;; 
	-1|--pass1)          nodraw=`toggle $nodraw`   ;; 
	-b|--batch)          batch=`toggle $batch`   ;; 
	-P|--proof)          proof=$2	      ; shift ;; 
	-M|--mc)             mc=`toggle $mc`   ;; 
	-g|--gdb)            gdb=`toggle $gdb`   ;; 
	-H|--hhd)            hhd=`toggle $hhd`   ;; 
	-O|--other)          comp=`toggle $comp`   ;; 
	-r|--rebin)          rebin=$2         ; shift ;; 
	-v|--vz-min)         vzmin=$2         ; shift ;; 
	-V|--vz-max)         vzmax=$2         ; shift ;; 
	-e|--energy)         cms=$2           ; shift ;;
	-S|--title)          tit="$2"         ; shift ;;
	-E|--eloss)          pass1=MakeELossFits.C ; 
	                     output=forward_eloss.root ;
			     nodraw=1 ;;
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
if test $nodraw -lt 1 ; then 
    opts1="-q" 
fi
if test $batch -gt 0 ; then 
    opts="-b -q $opts"
    redir="2>&1 | tee ${base}.log"
    echo "redir=$redir"
fi 
if test $noanal -lt 1 ; then 
    rm -f AnalysisResult.root AliAODs.root
    rm -f fmdana.png

    if test $gdb -gt 0 ; then 
	export PROOF_WRAPPERCMD="gdb -batch -x $ALICE_ROOT/PWG2/FORWARD/analysis2/gdb_cmds --args"
    fi
    echo "Running aliroot ${opts} ${opts1} ${ana}/${pass1}\(\".\",$nev,$proof,$mc\) $redir"
    if test $batch -gt 0 ; then 
	aliroot $opts $opts1 ${ana}/${pass1}\(\".\",$nev,$proof,$mc\) 2>&1 | tee ${base}.log
    else 
	aliroot $opts $opts1 ${ana}/${pass1}\(\".\",$nev,$proof,$mc\)
    fi
    fail=$?
    rm -f event_stat.root \
	EventStat_temp.root \
	outputs_valid \
	`printf %09d.stat $nev` 
    if  test $fail -gt 0  ; then 
        echo "Return value $fail not 0" ; exit $fail 
    fi
    if test ! -f ${output} ; then 
	echo "$output not made" ; exit 1; 
    fi
    if test ! -f AliAODs.root ; then 
	echo "No AOD creates" ; exit 1;
    fi
    echo "Analysis done"
fi

if test $nodraw -lt 1 ; then
    rm -f result.root 
    tit=`echo $tit | tr ' ' '@'` 
    echo "Running aliroot ${opts} ${opts1} ${ana}/Pass2.C\(\".\",\"$type\",$cms,$vzmin,$vzmax,$rebin,$hhd,$comp\)"
    aliroot ${opts} ${ana}/Pass2.C\(\".\",\"$type\",$cms,$vzmin,$vzmax,$rebin,\"$tit\",$hhd,$comp\)
fi


#
# EOF
#
