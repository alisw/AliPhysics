#!/bin/bash

# --------------------------------------------------------------------
uid=`id -u`
genv_file=/tmp/gclient_env_${uid}

if test ! -f ${genv_file} ; then 
    echo "No such file: ${genv_file}, please do alien-token-init" >/dev/stderr
    exit 1
fi
alien-token-info | grep -q "Token is still valid"
if test $? -ne 0 ; then 
    echo "Token not valid, please re-new" > /dev/stderr 
    exit 1
fi

# --------------------------------------------------------------------
verb=0

mess()
{
    if test $verb -lt 1 ; then return ; fi 
    echo $@
}

# --------------------------------------------------------------------
usage()
{
    cat <<EOF
Usage: $0 -p PRODUCTION [OPTIONS]

Options: 
	-h,--help		      This help 
	-v,--verbose                  Be verbose
	-p,--production  IDENTIFIER   Production identifier [$prod]
	-y,--year        YEAR         Year of production [$year]
	-P,--pass 	 NUMBER	      Reconstruction pass number [$pass]
	-d,--destination DIRECTORY    Directory to store result in [$dest] 
	-r,--run         NUMBER       Run number [$runn]
	-q,--qa          NUMBER       QA number 
	-a,--archives		      Get ZIP archives 
	-n,--no-action		      Run dry (do not copy files)

If run number is set, get the parts of next-to-last merge of that run only.

EOF
}

# --------------------------------------------------------------------
check_file()
{
    inp=$1 
    scr=`mktemp -u Check_XXXXXXXX` 
    cat <<EOF > ${scr}.C  
void ${scr}() 
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  // gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWG1");
  gSystem->Load("libPWG3base");
  TFile* file = TFile::Open("${inp}", "READ");
  if (!file) { 
    Error("${scr}", "No such file ${inp}");
    exit(1);
    return false;
  }
  TObject* forward1 = file->Get("Forward");
  if (!forward1) {
    Error("${scr}", "No Forward object found in ${inp}");
    exit(2);
    ret = false;
  } 
  TObject* forward2 = file->Get("ForwardResults");
  if (!forward2) {
    Error("${scr}", "No ForwardResults object found in ${inp}");
    exit(4);
  } 
  exit(0);
}
EOF
    # cat ${scr}.C 
    mess -n "aliroot -l -b -q ${scr}.C "
    aliroot -l -b -q ${scr}.C > /dev/null 2>&1 
    ret=$? 
    mess "-> $ret (0: good, 1: no file, 2: no Forward, 4: no ForwardResults"
    rm -f ${scr}.C 
    return $ret
}

# --------------------------------------------------------------------
year=
prod=
pass=1
qual=
dest=.
runn=0
qano=0
noac=0
arch=0

while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)        usage   ; exit 0 ;; 
	-v|--verbose)	  let verb=$verb+1  ;; 
	-p|--production)  prod=$2 ; shift ;; 
	-P|--pass)        pass=$2 ; shift ;; 
	-Q|--prepass)     qual=$2 ; shift ;; 
	-y|--year)        year=$2 ; shift ;; 
	-d|--destination) dest=$2 ; shift ;; 
	-r|--run)         runn=$2 ; shift ;; 
	-q|--qa)          qano=$2 ; shift ;; 
	-n|--no-action)   noac=1  ;; 
	-a|--archives)    arch=1  ;; 
	*) echo "$0: Unknown option $1" > /dev/stderr ; exit 2 ;; 
    esac
    shift
done

# --------------------------------------------------------------------
if test "x$prod" = "x" ; then 
    echo "No production identifier given" > /dev/stderr 
    exit 2
fi

if test "x$year" = "x" ; then 
    year=`echo $prod | sed -e 's/LHC\(..\).*/\1/'` 
    if test "x$year" = "x" ; then 
	echo "Couldn't get year from production identifier $prod" > /dev/stderr
	exit 2
    fi
fi

redir="/dev/null"
if test $verb -gt 1 ; then redir1=/dev/stderr ; fi

# --------------------------------------------------------------------
path=/alice/data/20${year}/${prod}/
store=${dest}/${prod}/${qual}pass${pass}
search="ESDs/${qual}pass${pass}/"
if test $runn -gt 0 ; then 
    path=`printf "${path}%09d/ESDs/${qual}pass${pass}/" $runn` 
    store=`printf "${store}/%09d" $runn` 
    search=
fi
if test $qano -gt 0 ; then 
    path=`printf "%sQA%02d/" $path $qano` 
fi
if test $arch -gt 0 ; then 
    search="${search}QA_archive.zip"
else
    search="${search}QAresults.root"
fi

cat <<EOF
Settings:

	Production:		$prod
	Year:			$year
	Pass:			$pass
	Pass qualifier:		$qual
	Path:                   $path
	Destination:		$dest
	Store:			$store
	Run number:             $runn
	Search string:		$search
EOF
    
# --------------------------------------------------------------------
mkdir -p ${store}
mess "alien_find ${path} ${search}"
files=`alien_find ${path} ${search} | grep -v "files found" 2> ${redir}` 
j=0
runs=
for i in $files ; do 
    b=`echo $i | sed -e "s,${path},,"` 
    d=
    if test $runn -gt 0 ; then 
	r=`echo $b | sed -e "s,[0-9]*${runn}\([0-9.]*\)/.*,\1," | tr '.' '_'`
	if test $arch -lt 1 ; then 
	    o=${store}/QAresults_${r}.zip
	else
	    d=${store}/${r}
	    mkdir -p $d
	    o=${d}/QAarchive.zip 
	fi
    else 
	r=`echo $b | sed -e "s,/.*,,"` 
	o=${store}/QAresults_${r}.root
    fi
    runs="$runs $r" 
    mess "$i -> $o"
    if test $noac -lt 1 && test ! -f $o ; then 
	mess "alien_cp alien:${i} file:${o}"
	alien_cp alien:${i} file:${o} > ${redir} 2>&1 
    fi
    if test $noac -lt 1 && test $arch -lt 1  ; then check_file $o ; fi
    if test $noac -lt 1 && test $arch -gt 0 ; then 
	(cd ${d} && unzip QAarchive.zip) > $redir 2>&1 
    fi
    let j=$j+1
done 
mess "Got a total of $j files for runs $runs"

# --------------------------------------------------------------------
#
# EOF
#

    