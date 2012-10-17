#!/bin/sh

type=
proto=
host=
file=
opts=
anchor="esdTree"
dbg=
ex=
par=0
dir=/data/alice/data/ppb/LHC12g/pass1/188359/

usage()
{
cat<<EOF
Usage: $0 -t TYPE [OPTIONS]

Options:
	-t,--type  TYPE       Set type of job [$type]
	-v,--verbose LEVEL    Set verbosity
	-d,--debug            Run in debugger [$dbg]
	-b,--batch            Run batch mode [$ex]
	-p,--par	      Run with PAR files [$par]
        -h,--help             Show this help

TYPE is one of 

	local        Run over local files 
	lite         Run in Proof-lite
	hehi         Run on hehi Proof farm
	caf          Run on CAF 
	caf plugin   Run on CAF using plugin 
	grid         Run in AliEn - Real data
	grid hijing  Run in AliEn - Hijing simulation
	grid dpmjet  Run in AliEn - DPMJet simulation
EOF
}

while test $# -gt 0 ; do 
    case $1 in 
	-t|--type) type=$2 ; shift ;; 
	-d|--debug) dbg="gdb --args" ;; 
	-b|--batch) ex="$ex --batch" ;; 
	-v|--verbose) ex="$ex --verbose=$2" ; shift ;; 
	-p|--par)   par=1 ;;
	-h|--help)  usage ; exit 0 ;;
    esac
    shift 
done
case $type in 
    local|lite) 
	proto=$type
	file="$dir"
	if test "x$type" = "xlite" ; then opts="workers=10" ; fi
	;;
    hehi*)
	proto=proof
	host="hehi00.nbi.dk";
	file="/default/cholm/LHC12g_pass1_uncalibrated_ESD_188359_partial";
	opts="mode=default&dsname"
	if test $par -gt 0 ; then opts="$opts&par=tasks" ; fi
	;;
    caf*)
	proto=proof
	host="alice-caf.cern.ch";
	file="/alice/data/LHC12g_000188359_ESDs_p1_uncalibrated";
	opts="mode=default&workers=60&aliroot=v5-03-68-AN&dsname"
	if test $par -gt 0 ; then opts="$opts&par=tasks" ; fi
	case $type in 
	    *plugin*) opts="${opts}&plugin" ;;
	esac
	;;
    grid*)
	proto=alien
	opts="run=188359-188359"
	case $type in 
	    *hijing*)
		file="/alice/sim/2012/LHC12g1/";
		opts="$opts&pattern=*/AliESDs.root&mc"
		;;
	    *dpmjet*)
		file="/alice/sim/2012/LHC12g4a/";
		opts="$opts&pattern=*/AliESDs.root&mc"
		;;
	    *)
		file="/alice/data/2012/LHC12g";
		opts="$opts&pattern=ESDs/pass1_uncalibrated/*/AliESDs.root"
		;;
	esac
	if test $par -gt 0 ; then opts="$opts&par=tasks" ; fi
	;;
    help*)
	;;
esac

echo "Building runTrain ..."
bld="g++ -g `root-config --cflags --glibs` \
    -lVMC -lGeom -lMinuit -lXMLIO -lTree -lTreePlayer \
    -I$ALICE_ROOT/include -L$ALICE_ROOT/lib/tgt_${ALICE_TARGET} \
    -lSTEERBase -lESD -lAOD -lANALYSIS -lOADB -lANALYSISalice \
    trainMain.cxx -o runTrain"
echo $bld
$bld || exit 1

echo "Cleaning previous run ($name) ..."
name=test_`echo $type | tr ' \t&/' '_'`
rm -rf $name

echo "Running the job ..."
cmd="./runTrain --class=MakeAODTrain --name=$name \
    --url=${proto}://${host}/${file}?${opts}#${anchor} \
    --sys=3 --snn=5023 --field=-5  --cent \
    $ex $@"
echo "$dbg $cmd"
$dbg $cmd

echo "done"
