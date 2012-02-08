#!/bin/sh 

nev=10000
noanal=0
nodraw=0
rebin=1
vzmin=-10
vzmax=10
type=
cms=900

usage()
{
cat<<EOF
Usage: $0 [OPTIONS]

Run analysis and visualize it using HHD analysis code 

Options:
	-h,--help		This help                  
	-N,--no-analysis	Do not analyse, just draw  ($noanal)
	-D,--no-draw		Do not draw, just analysis ($nodraw)
	-n,--events N		Number of events           ($nev)
	-r,--rebin N		Rebin by N                 ($rebin)
	-v,--vz-min CM          Minimum value of vz        ($vzmin)
	-V,--vz-max CM          Maximum value of vz        ($vzmax)
	-t,--trigger TYPE       Select trigger TYPE        ($type)
	-e,--energy CMS         Center of mass energy      ($cms)

TYPE is one of 
 
  INEL	
  INEL>0
  NSD 
EOF
}

while test $# -gt 0 ; do
    case $1 in 
	-h|--help)           usage            ; exit 0;; 
	-n|--events)         nev=$2           ; shift ;; 
	-N|--no-analysis)    noanal=$((($noanal+1)%2))   ;; 
	-D|--no-draw)        nodraw=$((($nodraw+1)%2))   ;; 
	-r|--rebin)          rebin=$2         ; shift ;; 
	-v|--vz-min)         vzmin=$2         ; shift ;; 
	-V|--vz-max)         vzmax=$2         ; shift ;; 
	-e|--energy)         cms=$2           ; shift ;;
	-t|--type)           type="$2"        ; shift ;;
	*) echo "$0: Unknown option '$1'" >> /dev/stderr ; exit 1 ;;
    esac
    shift
done 

base=`printf hhd_%07d $nev` 

keep_file()
{
    file=$1
    if test -f $file ; then 
	mv $file ${file}.keep 
    fi
}
restore_file()
{
    file=$1
    if test -f ${file}.keep ; then 
	mv ${file}.keep ${file}
    fi
}

if test $noanal -lt 1 ; then 
    keep_file AnalysisResult.root
    keep_file AliAODs.root 

    aliroot -l -b -q -x RunManager.C\(\".\",$nev\) 2>&1 | tee $base.log 

    if test ! -f AnalysisResults.root ; then 
	echo "Analysis failed" 
	restore_file AnalysisResult.root
	restore_file AliAODs.root 
	exit 1
    fi
    echo "Analysis done"
    rm -f event_stat.root EventStat_temp.root outputs_valid
    mv AnalysisResults.root ${base}_hists.root
    if test -f AliAODs.root ; then 
	mv AliAODs.root     ${base}_aods.root
    fi
    restore_file AnalysisResult.root
    restore_file AliAODs.root 
fi


if test $nodraw -lt 1 ; then
    rm -f fmd_dNdeta_mult.root

    samp=0
    tt=`echo $type | tr '[a-z]' '[A-Z]'` 
    case "x$tt" in 
	"xINEL")   n="inel";    s=""    ; samp=0 ;; 
	"xINEL>0") n="inelgt0"; s="_NSD"; samp=4 ; echo "Using NSD for $tt";; 
	"xNSD")    n="nsd";     s="_NSD"; samp=4 ;; 
	*) echo "Unknown type $tt" > /dev/stderr ; exit 1 ;; 
    esac
	    
    aliroot -l -b -q -x \
	$ALICE_ROOT/PWGLF/FORWARD/analysis/drawdNdeta.C\(\"${base}_hists.root\",$samp,$rebin,$vzmin,$vzmax,1,$cms\)
    tmin=`printf %+03d $vzmin | tr '+-' 'pm'` 
    tmax=`printf %+03d $vzmax | tr '+-' 'pm'` 
    out=`printf "hhd_%04dGeV_%s-%scm_rb%02d_%s.root" $cms $tmin $tmax $rebin $n`
    mv fmd_dNdeta_mult${s}.root $out

    mv fmdana.png `basename $out .root`.png 
fi
#
# EOF
#
