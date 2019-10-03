#!/bin/bash
#
# variables to be set 
#
#  base:         Location of base scripts
#  datadir:      Top of data directories 
#  noact:        Whether to do anything at all 
#

# --- Test no-act flag -----------------------------------------------
isNoAction()
{
    if test "x$noact" != x && test $noact -gt 0 ; then 
	true 
    else
	false
    fi 
}

# --- Write configuration file ---------------------------------------
writeConfig()
{
    if isNoAction ; then return ; fi 

    local sflow=$1 ; shift 
    local sfhigh=$1; shift 
    local dclow=$1 ; shift
    local three=$1 ; shift
    sed -e "s/@sflow@/${sflow}/" \
	-e "s/@sfhigh@/${sfhigh}/" \
	-e "s/@dclow@/${dclow}/" \
	-e "s/@three@/${three}/" \
	< ${base}/ForwardAODConfig.C.in > ForwardAODConfig.C 
}

# --- Get the URL ----------------------------------------------------
# mc=0
getUrl()
{
    local run=$1 ; shift 
    local nwrk=8
    local dir=${datadir}/${run}/index.root
    local opts="workers=${nwrk}"
    if test x$mc != x  ; then 
	opts="${opts}&mc"
    fi
    local url="lite://$dir?${opts}#esdTree"
    echo "$url"
}

# --- Get cut string -------------------------------------------------
getCut()
{
    local met=`echo $1 | tr '[A-Z]' '[a-z]'` ; shift
    local cut=`echo $1 | tr '[A-Z]' '[a-z]'` ; shift
    local mnm="kMPVFraction"
    local rcut=$cut 
    local same=1

    # --- Get what we need to do -------------------------------------
    case $met in 
	mpv)  mnm=kMPVFraction ;;
	xi)   mnm=kLandauWidth ;;
	sig)  mnm=kLandauSigmaWidth ;;
	fix)  mnm=kFixed ;;
	fit)  mnm=kFitRange ;;
	prob) mnm=kProbability ; rcut="1e-${cut}";;
	*)   echo "Unknown method type: $met" ; exit 1;;
    esac
    ret="${mnm},${rcut}"
    if test $same -gt 0 ; then 
	ret="$ret,${rcut},${rcut},${rcut},${rcut}"
    fi
    echo "$ret"
}
log=log

# --- Print outs -----------------------------------------------------
debug=2
msg()
{
    local lvl=$1 ; shift 
    if test $lvl -le $debug ; then 
	printf "%s\n" "$*"
    fi
}

# --- Run one job ----------------------------------------------------
iter=1
niter=1
toDo=0
runOne()
{
    local run=`echo $1 | tr '[A-Z]' '[a-z]'`  ; shift 
    local thr=`echo $1 | tr '[A-Z]' '[a-z]'`  ; shift 
    local slmet=`echo $1 | tr '[A-Z]' '[a-z]'`; shift
    local slcut=`echo $1 | tr '[A-Z]' '[a-z]'`; shift
    local shmet=`echo $1 | tr '[A-Z]' '[a-z]'`; shift
    local shcut=`echo $1 | tr '[A-Z]' '[a-z]'`; shift
    local dcmet=`echo $1 | tr '[A-Z]' '[a-z]'`; shift
    local dccut=`echo $1 | tr '[A-Z]' '[a-z]'`; shift

    case x$thr in 
	x)             three=true  ; nstr=3strip ;;
	x|xfalse|xno*) three=false ; nstr=2strip ;;
	*)             three=true  ; nstr=3strip ;;
    esac

    # --- Set up some variables --------------------------------------
    local url=`getUrl $run`
    local class="MakeAODTrain"
    local dcWhat=`printf "%s%05.2f" $dcmet $dccut | tr '.' 'd'` 
    local slWhat=`printf "%s%05.2f" $slmet $slcut | tr '.' 'd'` 
    local shWhat=`printf "%s%05.2f" $shmet $shcut | tr '.' 'd'` 
    local what="${nstr}_sl${slWhat}_sh${shWhat}_dc${dcWhat}"
    local name="${run}_aod_${what}"

    # --- Check if we already have the output ------------------------
    exists=
    if test -d ${name} && test -s ${name}/AliAOD.root; then 
	# echo "${name} already exists" 
	exists=" done"
    else
	let toDo=$toDo+1
    fi

    local prog=`printf "%3d/%3d" $iter $niter`
    msg 1 "=         ${prog} ${name}${exists}" 
    let iter=$iter+1
    if test "X${exists}" != "X" ; then 
	return
    fi

    # --- Write configuration ----------------------------------------
    dcSub=`getCut "$dcmet" "$dccut"` 
    slSub=`getCut "$slmet" "$slcut"` 
    shSub=`getCut "$shmet" "$shcut"` 
    # echo "=         Write config w/ sl=$slSub, sh=$shSub, dc=$dcSub and $nstr"
    writeConfig \
	"AliFMDMultCuts::$slSub" \
	"AliFMDMultCuts::$shSub" \
	"AliFMDMultCuts::$dcSub" \
	"$three"

    # --- Our options ------------------------------------------------
    opts=(--class=$class \
	--name=$name \
	--cent \
	--central-config="CentralAODConfig.C" \
	--corr="." \
	--date="none" \
	--forward-config="ForwardAODConfig.C" \
	--type="ESD" \
	--url="$url")
	
    # --- Remove old directory ---------------------------------------
    if ! isNoAction ; then 
	rm -rf ${name}
	mkdir -p ${log}
    fi

    # --- Run the job ------------------------------------------------
    msg 3 "=         Running runTrain ${opts[@]} $@ - see log/$name.log"
    if ! isNoAction ; then 
	runTrain "${opts[@]}" $@ > log/${name}.log 2>&1 
	# --- Wait for things to settle ------------------------------
	sleep 1
    fi

    
    # --- Check if we succeeded, and summarize -----------------------
    if ! isNoAction ; then 
	if test -f ${name}/forward.root && \
	    test ! -f ${name}/forward.pdf ; then 
	    msg 3 "=         Post processing in ${name}"
	    (cd ${name} && ./post.sh) >> log/${name}.log 2>&1 
	fi
    else
	msg 3 "=         Post processing in ${name}"
    fi

    # --- Check if we got an AOD out, and run dN/deta analysis -------
    dname=`echo $name | sed 's/aod/dndeta/'` 
    if ! isNoAction ; then 
	if test -f ${name}/dndeta.sh && test -s ${name}/AliAOD.root ; then 
	    rm -rf ${dname}
	    msg 3 "=         Running ./${name}/dndeta.sh  - see log/$dname.log"
	    ./${name}/dndeta.sh --batch $@ > log/$dname.log 2>&1 
	    
	    # --- Wait for things to settle --------------------------
	    sleep 1
	    
	    # --- Check if we succeeded, and summarize ---------------
	    if test -f ${dname}/forward_dndeta.root && \
		test ! -f ${dname}/forward_dndeta.pdf ; then 
		msg 3 "=         Post processing in ${dname}"
		(cd ${dname} && ./post.sh) >> log/$dname.log 2>&1 
		fi
	    fi
    else 
	msg 3 "=         Running ./${name}/dndeta.sh  - see $dname.log"
	msg 3 "=         Post processing in ${dname}"
    fi
}
# --- Extract the method from string ---------------------------------
extractMethod()
{
    echo "$1" | cut -f1 -d=
}
# --- Extract the parameters from string -----------------------------
extractPars()
{
    echo "$1" | cut -f2 -d= | tr ',' ' '
}
# --- Extract the cut ------------------------------------------------
extractWhich()
{
    echo "$1" | cut -f1 -d:
}
# --- Count the number of cuts to loop over --------------------------
cntCuts()
{
    local c="$1" ; shift 
    local w=`extractWhich $c`
    local l=`echo "$c" | cut -f2 -d: | tr ';' ' '`
    local n=0

    for i in $l ; do 
	local m=`extractMethod $i`
	local p=`extractPars $i` 
	for j in $p ; do 
	    let n=$n+1
	done
    done
    echo $n
}

# --- Loop over cuts -------------------------------------------------
loopCuts()
{
    local r="$1" ; shift 
    local t="$1" ; shift 
    local q="$1" ; shift 
    local a="$1" ; shift 
    local c="$1" ; shift 
    local e=0
    if test "X$1" == "X--" ; then 
	e=1 
	shift
    fi
    # echo "r=$r t=$t q=$q a=$a c=$c @=$@"
    local w=`extractWhich $c`
    local l=`echo "$c" | cut -f2 -d: | tr ';' ' '`
    # echo "Which: $w Cuts: $l"
    
    
    for i in $l ; do 
	local m=`extractMethod $i`
	local p=`extractPars $i` 
	msg 2 "${q}$w:${m} ($p)"
	for j in $p ; do 
	    # echo "j=$j e=$e"
	    if test $e -gt 0 ; then 
		msg 2 "${q} W/cut $w:$m:$j ($a)"
		runOne $r $t $a $m $j $@ 
		msg 2 "${q} done w/cut $w:$m:$j"
	    else
		msg 2 "${q}Loop w/cut $w:$m:$j ($a, $@)"
		loopCuts $r $t "${q}  " "${a} ${m} ${j}" "$@"
		msg 2 "${q}done loop w/cut $w:$m:$j"
	    fi
	done
    done
}
    
# --- Extract to script ----------------------------------------------
writeScriptCut() 
{ 
    local spec=$1 ; shift 
    local which=`extractWhich $spec | tr '[a-z]' '[A-Z]'`
    local l=`echo "$spec" | cut -f2 -d: | tr ';' ' '` 
    # echo "which=$which spec=$spec l=$l" > /dev/stderr
    for i in $l ; do 
	local m=`extractMethod $i` 
	local p=`extractPars   $i` 
	# echo " i=$i m=$m p=$p" > /dev/stderr 
	echo "  t.Add${which}Cut(\"$m\",\"$p\");" 
    done
}

# --- The full loop --------------------------------------------------
fullLoop()
{
    local runs="$1"       ; shift 
    local strcuts="$1"    ; shift 
    local slcuts="sl:$1"  ; shift 
    local shcuts="sh:$1"  ; shift 
    local dccuts="dc:$1"  ; shift 
    
    local nrn=`echo $runs | wc -w` 
    local nst=`echo $strcuts | wc -w` 
    local nsl=`cntCuts $slcuts`
    local nsh=`cntCuts $shcuts`
    local ndc=`cntCuts $dccuts`
    local nto=$((${nrn}*${nst}*${nsl}*${nsh}*${ndc}))
    niter=$nto

    cat <<EOF
runs:      $nrn ${runs}
strcuts:   $nst ${strcuts}
slcuts:    $nsl ${slcuts}
shcuts:    $nsh ${shcuts}
dccuts:    $ndc ${dccuts}
rest:      $@ 

Number of iterations: ${nrn} * ${nst} * ${nsh} * ${ndc} = ${nto}
EOF
    cat <<EOF > Trending.C
void Trending() { 
  const char* fwd = "\$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  gROOT->SetMacroPath(Form("%s:%s:%s/scripts:\$(ANA_SRC)/scan", 
			   gROOT->GetMacroPath(), 
			   fwd, fwd));
  gSystem->AddIncludePath(Form("-I%s", fwd));
  gSystem->AddIncludePath(Form("-I%s/scripts", fwd));

  gROOT->LoadMacro("Trend.C++g");

  Trend t;
EOF
    writeScriptCut "${slcuts}" >> Trending.C
    writeScriptCut "${shcuts}" >> Trending.C
    writeScriptCut "${dccuts}" >> Trending.C
    for r in $runs ; do
	echo "  t.AddRun(\"$r\");" >> Trending.C
    done 
    cat <<EOF >> Trending.C
  t.AddCentrality( 0, 5);
  t.AddCentrality( 5,10);
  t.AddCentrality(10,20);
  t.AddCentrality(20,30);
  t.SetOrder("sl sh dc");

  TString out("scanProbXiSigma.root");
  t.Run(out);
}
EOF

    # --- Loop over everyting --------------------------------------------
    for r in $runs ; do                      # Loop over runs 
	msg 2 "= Processing run $r" 
	for t in ${strcuts} ; do             # Loop over 3-strip merging 
	    msg 2 "=  Allow 3-particle merging $t" 
	    
	    loopCuts "$r" "$t" "=   " "" "$slcuts" "$shcuts" "$dccuts" -- $@ 
	done
	msg 2 "= Done with run $r"
    done                                         # Loop over runs 
    echo "Did $toDo iterations out of $niter"
}

# EOF
