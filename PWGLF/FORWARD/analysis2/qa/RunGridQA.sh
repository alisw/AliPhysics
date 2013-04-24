#!/bin/bash

# --- Usage ----------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
	-p,--production PERIOD      LHC Period identifier [$prod]
	-P,--pass       PASS        Production pass identifier [$prodpass]
	-n,--part       PART        Part identifier [$part]
	-r,--runs       RUNS        List of runs or file [$runs]
	-c,--clean	  	    Clean previous job output [$clean]
	-v,--verbose		    Increase verbosity [$verb]
EOF

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

# --- Diagnostics output ---------------------------------------------
verb=0
mess()
{
    if test $1 -gt $verb ; then return ; fi 
    shift
    echo $*
}

# --- Parse production information -----------------------------------
parse_prod()
{
    prodyear=`echo $prodfull | sed 's/LHC\(..\).*/\1/'` 
    prodletter=`echo $prodfull | sed "s/LHC${prodyear}\(.\).*/\1/"` 
    prodpost=`echo $prodfull | sed "s/LHC${prodyear}${prodletter}//"` 
}

parse_pass()
{
    passno=`echo $passfull | sed 's/.*pass\([0-9]*\).*/\1/'`  
    passpost=`echo $passfull | sed -n "s/.*pass${passno}_//p"` 
    passpre=`echo $passfull | sed -n "s/pass.*//p"` 
}

# --- Append path element --------------------------------------------
append_to_path()
{
    local tmp=$1 ; shift 
    local add=$1
    if test "x$tmp" != "x" ; then tmp="${tmp}/" ; fi 
    echo ${tmp}${add}
}

# --- Get the path and search pattern --------------------------------
path=
search=
setup_input()
{
    local datd=data
    local esdd=ESDs/
    case x$prodpost in 
	x_*) ;; 
	x) ;; 
	*)  mess 3 "Assuming simulation output"
	    datd=sim 
	    esdd= 
	    ;; 
    esac
    
    local paid=
    if echo "$passno" | grep -q -E '^[0-9]*[.]?[0-9]*$' ; then 
	if test "x$passfull" != "x" && test $passno -gt 0 ; then 
	    paid=pass${passno}
	fi
    else
	paid=${passfull}
	passpre=
	post=
    fi
    local post=${passpost}
    case x$post in 
	x_*) ;; 
	x) ;; 
	*) post="_${post}" ;; 
    esac

    # Assume official productions 
    path=/alice/${datd}/${year}/${prodfull}/
    search="${esdd}${passpre}${paid}${post}"
    search=`append_to_path "$search" "*/AliESDs.root"` 
}

# --- Setup the runs -------------------------------------------------
unique=
setup_runs()
{
    if test -f $runs ; then 
	unique=`basename $runs .txt` 
	unique=`basename $unique .list` 
	unique=`basename $unique .runs` 
	return
    fi
    local l=`echo "$runs" | tr ',+:.' ' '` 
    local first=
    local last=
    for i in $l ; do 
	if test "x$first" = "x" ; then first=$i ; fi 
	last=$l 
    done
    unique=${first}_${last} 
    runs=`echo "$l" | tr ' ' ','` 
}

# --- Clean previous attempt -----------------------------------------
clean_previous()
{
    local alienu=`alien_whoami | tr -d ' '` 
    local alienl=`echo ${alienu} | cut -c1` 
    local aliend="/alice/cern.ch/user/${alienl}/${alienu}/${name}"
    mess 1 "Removing local directory ${name}"
    if test $noact -lt 1 ; then rm -rf ${name} ; fi
    mess 1 "Removing alien directory ${aliend}"
    if test $noact -lt 1 ; then alien_rmdir ${aliend} ; fi
}

# --- Parse command line options -------------------------------------
redir=/dev/null
part=
clean=0
noact=0
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)        usage ; exit 0 ;; 
	-v|--verbose)     let verb=$verb+1   ;; 
	-l|--log-file)    redir=             ;; 
	-p|--production)  prodfull=$2; shift; parse_prod ; year=20${prodyear} ;;
	-P|--pass)        passfull=$2; shift; parse_pass ;;
	-n|--part)        part=$2 ; shift ;; 
	-r|--runs)        runs=$2		; shift ;;
	-c|--clean)       clean=1 ; shift ;; 
	-s|--noact)       noact=1 ; shift ;;
	*) echo "$0: Unknown argument: $1" > /dev/stderr ; exit 1 ;; 
    esac
    shift
done

# --- Initial checks -------------------------------------------------
if test "x$prodfull" = "x" ; then 
    echo "No production specified" > /dev/stderr 
    exit 1
fi
if test "x$passfull" = "x" ; then 
    echo "No pass specified" > /dev/stderr 
    exit 1
fi
if test "x$runs" = "x" ; then 
    echo "No runs specified" > /dev/stderr 
    exit 1
fi

check_token 
setup_input
setup_runs
name="QA_${prodfull}_${passfull}"
if test "x$part" != "x" ; then 
    if echo "$part" | grep -q -E '^[0-9]*[.]?[0-9]*$' ; then 
	part="part${part}"
    fi
    name="${name}_${part}"
fi

if test $clean -gt 0 ; then clean_previous ; fi 
    

# --- Some friendly information --------------------------------------
cat <<EOF
	Year:			${year}
	Production:		${prodfull} 
	  Year:			${prodyear}
	  Letter:		${prodletter}
	  Suffix:		${prodpost}
	Pass:			${passfull}
	  Number:		${passno}
	  Prefix:		${passpre}
	  Postfix:		${passpost}
	Lock file:		${lock}
	Log:                    ${redir}
	Runs:                   ${runs}
	Path:                   ${path}
	Search:                 ${search}
	Part:                   ${part}
	Name:                   ${name}
EOF

# --- Do the actual running ------------------------------------------
url="alien:///alice/data/${year}/${prodfull}\?run=${runs}\&pattern=${search}#esdTree"
train=MakeQATrain
opts=(--class=${train} \
    --name=${name} \
    --cent \
    --url="${url}")
mess 1 "Will do runTrain \"${opts[@]}\" $@" 
if test $noact -gt  0 ; then 
    sleep 1 
else 
    runTrain "${opts[@]}" $@
fi

mess 1 "Will monitor in directory $name" 
if test $noact -gt 0 ; then 
    sleep 1 
else 
    (cd ${name} && aliroot -b -x -q -l Watch.C)
fi

#
# EOF
#


