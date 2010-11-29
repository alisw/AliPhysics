#!/bin/bash 

nev=10000
noanal=0
nodraw=0
rebin=1
vzmin=-10
vzmax=10
ncutbins=1
correctioncut=0.1
batch=0
gdb=0
proof=0
type=
cms=900
hhd=1
comp=1

usage()
{
cat<<EOF
Usage: $0 [OPTIONS]

Do Pass1 and Pass2 on ESD files in current directory.  

Options:
	-h,--help		This help                  
	-n,--events N		Number of events           ($nev)
	-N,--no-analysis	Do not analyse, just draw  ($noanal)
	-D,--no-draw		Do not draw, just analysis ($nodraw)
	-r,--rebin N		Rebin by N                 ($rebin)
	-c,--n-cut-bins N       Number of cut bins         ($ncutbins)
	-C,--correction-cut V   Cut on secondary corr,     ($correctioncut)
	-v,--vz-min CM          Minimum value of vz        ($vzmin)
	-V,--vz-max CM          Maximum value of vz        ($vzmax)
	-t,--trigger TYPE       Select trigger TYPE        ($type)
	-e,--energy CMS         Center of mass energy      ($cms)
	-b,--batch              Do batch processing        ($batch)
	-p,--proof		Run in PROOF(Lite) mode    ($proof)
	-g,--gdb		Run in GDB mode    	   ($gdb)
	-H,--hhd		Do comparison to HHD	   ($hhd)
	-O,--other		Do comparison to other	   ($comp)

TYPE is a comma or space separated list of 
 
  INEL	      Inelastic triggers (V0A|V0C|SPD)
  INEL>0      As above + N_ch > 0 in -0.5<eta<+0.5
  NSD         Non-single diffractive ((VOA&VOC)|N_ch > 5 -1.9<eta<+1.9)

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
	-N|--no-analysis)    noanal=`toggle $noanal`   ;; 
	-D|--no-draw)        nodraw=`toggle $nodraw`   ;; 
	-b|--batch)          batch=`toggle $batch`   ;; 
	-p|--proof)          proof=`toggle $proof`   ;; 
	-g|--gdb)            gdb=`toggle $gdb`   ;; 
	-H|--hhd)            hhd=`toggle $hhd`   ;; 
	-O|--other)          other=`toggle $other`   ;; 
	-r|--rebin)          rebin=$2         ; shift ;; 
	-v|--vz-min)         vzmin=$2         ; shift ;; 
	-V|--vz-max)         vzmax=$2         ; shift ;; 
	-e|--energy)         cms=$2           ; shift ;;
	-c|--n-cut-bins)     ncutbins=$2      ; shift ;; 
	-C|--correction-cut) correctioncut=$2 ; shift ;;
	-t|--type)           
	    if test "x$type" = "x" ; then type=$2 ; else type="$type|$2"; fi
	    shift ;;
	*) echo "$0: Unknown option '$1'" >> /dev/stderr ; exit 1 ;;
    esac
    shift
done 

base=`printf dndeta_%07d $nev`
opts="-l -x"
redir=
if test $batch -gt 0 ; then 
    opts="-l -b -q -x" 
    redir="2>&1 | tee ${base}.log"
fi 
if test $noanal -lt 1 ; then 
    rm -f AnalysisResult.root AliAODs.root
    rm -f fmdana.png

    if test $gdb -gt 0 ; then 
	export PROOF_WRAPPERCMD="gdb -batch -x $ALICE_ROOT/PWG2/FORWARD/analysis2/gdb_cmds --args"
    fi
    aliroot $opts $ALICE_ROOT/PWG2/FORWARD/analysis2/Pass1.C\(\".\",$nev,$ncutbins,$correctioncut,$proof\) $redir 
    rm -f event_stat.root EventStat_temp.root outputs_valid
    if test ! -f AnalysisResults.root || test ! -f AliAODs.root ; then 
	echo "Analysis failed" 
	exit 1
    fi
    echo "Analysis done"
fi

if test $nodraw -lt 1 ; then
    rm -f result.root 
    aliroot ${opts} $ALICE_ROOT/PWG2/FORWARD/analysis2/Pass2.C\(\"AliAODs.root\",\"$type\",$cms,$vzmin,$vzmax,$rebin,\"Run\ 118506,\ $nev\ Events,\ v_{z}#in[$vzmin,$vzmax],\ $type\",$hhd,$comp\)
fi


#
# EOF
#
