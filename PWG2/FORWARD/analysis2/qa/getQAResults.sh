#!/bin/bash
#
# @file   getQAResults.sh
# @author Christian Holm Christensen <cholm@nbi.dk>
# @date   Thu Nov 17 11:47:14 2011
# 
# @brief Retrieve trending.root/QAResults.root files from AliEn for a
# given producton as specified by the command line options 
# 
# @ingroup pwg2_forward_qa_scripts
#

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
	-P,--pass 	 NUMBER	      Reconstruction pass number [$pass]
	-Q,--pre-pass	 STRING	      Prefix to pass identifier [${prep}]
	-R,--post-pass	 STRING       Postfix to pass identifier [$post]
	-y,--year        YEAR         Year of production [$year]
	-d,--destination DIRECTORY    Directory to store result in [$dest] 
	-r,--run         NUMBER       Run number [$runn]
	-q,--qa          NUMBER       QA number 
	-f,--file	 NAME	      File to get [$file]
	-T,--trending	              Get trending.root file 
	-a,--archives		      Get ZIP archives 
	-n,--no-action		      Run dry (do not copy files)
	-i,--no-check		      Do not check files [$docheck]

If run number is set, get the parts of next-to-last merge of that run only.

EOF
}

# --------------------------------------------------------------------
docheck=1
check_file()
{
    if test $docheck -lt 1 ; then return 0; fi 
    inp=$1 
    scr=`mktemp -u Check_XXXXXXXX` 
    cat <<EOF > ${scr}.C  
void ${scr}() 
{
  int ret = 0;
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
  }
  TObject* forward1 = file->Get("Forward");
  if (!forward1) {
    Error("${scr}", "No Forward object found in ${inp}");
    ret |= 2;
  } 
  TObject* forward2 = file->Get("ForwardResults");
  if (!forward2) {
    Error("${scr}", "No ForwardResults object found in ${inp}");
    ret |= 4;
  } 
  exit(ret);
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
year=""
prod=""
pass=-1
prep=""
post=""
dest=.
runn=0
qano=0
noac=0
arch=0
file=QAresults.root 

while test $# -gt 0 ; do 
    case $1 in 
	-h|--help)        usage   ; exit 0 ;; 
	-v|--verbose)	  let verb=$verb+1  ;; 
	-p|--production)  prod=$2 ; shift ;; 
	-P|--pass)        pass=$2 ; shift ;; 
	-Q|--prepass)     prep=$2 ; shift ;; 
	-R|--postpass)    post=$2 ; shift ;; 
	-y|--year)        year=$2 ; shift ;; 
	-d|--destination) dest=$2 ; shift ;; 
	-r|--run)         runn=$2 ; shift ;; 
	-q|--qa)          qano=$2 ; shift ;; 
	-f|--file)        file=$2 ; shift ;; 
	-T|--trending)	  file=trending.root ; shift ;; 
	-a|--archives)    arch=1  ;; 
	-n|--no-action)   noac=1  ;; 
	-i|--no-check)    docheck=0 ;; 
	*) echo "$0: Unknown option $1" > /dev/stderr ; exit 2 ;; 
    esac
    shift
done

# --------------------------------------------------------------------
uid=`id -u`
genv_file=/tmp/gclient_env_${uid}

if test ! -f ${genv_file} ; then 
    echo "No such file: ${genv_file}, please do alien-token-init" >/dev/stderr
    exit 1
fi
. ${genv_file}
alien-token-info | grep -q "Token is still valid"
if test $? -ne 0 ; then 
    echo "Token not valid, please re-new" > /dev/stderr 
    exit 1
fi


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
lett=`echo $prod | sed -e "s/LHC${year}\(.\).*/\1/"`
suff=`echo $prod | sed -e "s/LHC${year}${lett}//"`

redir="/dev/null"
if test $verb -gt 1 ; then redir=/dev/stderr ; fi

# --------------------------------------------------------------------
if test "x$post" != "x" ; then 
    case $post in 
	_*) ;; 
	*) post="_${post}" ;; 
    esac
fi
if test ${pass} -ge 0 ; then 
    paid=pass${pass} 
fi
datd=data
esdd=ESDs/
if test "x${suff}" != "x" ; then 
    datd=sim
    esdd=
fi
path=/alice/${datd}/20${year}/${prod}/
store=${dest}/${prod}/${prep}${paid}${post}
search="${esdd}${prep}${paid}${post}"
if test $runn -gt 0 ; then 
    path=`printf "${path}%09d/ESDs/${prep}${paid}${pass}${post}/" $runn` 
    store=`printf "${store}/%09d" $runn` 
    search=
fi
if test $qano -gt 0 ; then 
    if test $runn -gt 0 ; then 
	path=`printf "%sQA%02d/" $path $qano` 
    else
	if test "x$search" != "x" ; then search=${search}/ ; fi
	search=${search}QA`printf %02d ${qano}` 
    fi
fi
if test $arch -gt 0 ; then 
    file=QA_archive.zip
fi
base=`basename $file .root`
if test "x$search" != "x" ; then search=${search}/ ; fi
search="${search}${file}"

# --------------------------------------------------------------------
cat <<EOF
Settings:

	Production:		$prod
	Year:			$year
	Letter:                 $lett
	Suffix:			$suff
	Pass:			$pass
	Pass prefix:		$prep
	Pass postfix:		$post
	File:			$file
	Path:                   $path
	Destination:		$dest
	Store:			$store
	Run number:             $runn
	Search string:		$search
	Verbosity:		$verb
	Redirection:		$redir
EOF
    
# --------------------------------------------------------------------
mkdir -p ${store}
mess "Getting list of files from AliEn - can take minutes - be patient"
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
	    o=${store}/${base}_${r}.root
	else
	    d=${store}/${r}
	    mkdir -p $d
	    o=${d}/${file}
	fi
    else 
	r=`echo $b | sed -e "s,/.*,,"` 
	o=${store}/${base}_${r}.root
    fi
    runs="$runs $r" 
    mess "$i -> $o"
    if test $noac -lt 1 && test ! -f $o ; then 
	mess "alien_cp alien:${i} file:${o}"
	alien_cp alien:${i} file:${o} > ${redir} 2>&1 
	if test -f $o ; then chmod g+rw ${o} ; fi
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

    
