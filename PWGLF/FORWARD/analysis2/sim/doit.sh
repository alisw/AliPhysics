#!/bin/bash
version=5
tag=
id=
run=
stage=0
upload=0 
jobs=2
events=1
gen=default
aliroot="v5-04-Rev-20"
root=""
geant=""
minmerge=30
noact=0

# --- Display help message -------------------------------------------
usage()
{
    cat <<EOF 
Usage: $0 [OPTIONS]

Options:
	-h,--help		This help
	-t|--tag       TAG	Job tag [****] ($tag)
	-i|--id        NAME	Name of production ($id)
	-R|--run       RUN_NO	Run number ($run)
	-e|--generator NAME	Event generator ($gen)
	-c|--copy		Copy files to AliEn
	-n|--jobs      JOBS	Set number of jobs[**] ($jobs)
	-m|--events    EVENTS	Set events/job[**] ($events)
	-s|--stage     STAGE	Set the stage[***] ($stage)
	-o|--output    DIR	Set base output directory ($aout)
	-d|--data      DIR      Set data directory ($adir)
	-b|--bin       DIR	Set base bin directory ($abin)
	-a|--aliroot   RELEASE	Set AliROOT release [*] ($aliroot)
	-r|--root      RELEASE	Set ROOT release [*] ($root)
	-g|--geant     RELEASE	Set GEANT3 release [*] ($geant)
	-f|--final     NUMBER 	Run final merging when down to this ($minmerge)

[*] Only make sense with option -c 
[**] Only make sense for stage 0
[***] If stage is set to 6, try to deduce the stage automatically
[****] TAG is a short hand for specific id and run
EOF
}

# --- Process a return value -----------------------------------------
log_msg() 
{
    local log=$1 ; shift
    echo -en "$@\e[0m ... "
    if test "x$log" != "x" ; then 
	echo "=== $@" >> $log 
    fi
}
# --- Make error -----------------------------------------------------
log_err() 
{
    local pre=$1 
    local post=$2 
    echo -e "\e[31;1mError\e[0m: ${pre} \e[35m${post}\e[0m" > /dev/stderr 
}

# --- Process a return value -----------------------------------------
log_end()
{
    local log=$1
    local ret=$2 
    local ext=$3
    local msg=""
    local fmt=""
    if test $ret -eq 0 ; then 
	msg="success"
	fmt="32"
    else 
	msg="failure"
	fmt="31"
    fi 
    echo -e "\e[${fmt}m${msg}${ext}\e[0m"
    if test "x$log" != "x" ; then 
	echo "=== $msg$ext" >> $log 
    fi 
}

# --- Copy a file to AliEn -------------------------------------------
copy()
{
    local file=$1 
    local dest=$2
    local del=$3
    local base=`basename $file`

    if test "x$del" != "x" ; then 
	log_msg cp.log "Removing \e[33malien:${dest}/${base}"
	if test $noact -lt 1 ; then 
	    alien_rm ${dest}/${base} >> cp.log 2>&1 
	fi 
	log_end cp.log $? " (ignore errors)"
    fi

    log_msg cp.log "Uploading \e[33m${file}\e[0m to \e[33m${dest}"
    if test $noact -lt 1 ; then 
	alien_cp -n file:${file} alien:${dest}/${base} >> cp.log 2>&1
    fi 
    log_end cp.log $? 
}

# --- Run merging jpb ------------------------------------------------
merge()
{
    local what=$1
    local stage=$2 
    local dir=$3
    local out=$4
    local tag=$5 
    local run=$6
    local tmpdir=`mktemp -d` 
    
    local top=${out}/${tag}/${run}
    local bse=${what}_Stage_${stage}.xml 
    local xml=${tmpdir}/${bse}
    local arc=${what}_archive.zip
    local jdl=Merge.jdl
    local ret=0

    rm -f cp.log 

    log_msg cp.log "Creating XML file \e[33m${xml}"
    if test $stage -eq 1 ; then 
	rm -f ${xml}
	alien_find -x ${top}/${bse} ${top} */${arc} > ${xml} 2>>cp.log
	ret=$? 
    else 
	let prev=$stage-1
	rm -f ${xml}
	alien_find -x ${top}/${bse} ${top}/${what}_Stage_${prev} */${arc} \
	    > ${xml} 2>>cp.log
	ret=$? 
    fi
    log_end cp.log $ret 
    if test $ret -ne 0 ; then 
	log_err "Make XML", "Failed to make XML collection $bse"
	exit 1
    fi 
    local n=`grep "<event name" ${xml} | wc -l 2>/dev/null` 
    if test $n -lt $minmerge ; then 
	old=$bse
	stage=5
	jdl=Final.jdl
	bse=${what}_Stage_${stage}.xml 
	tmp=${tmpdir}/${bse}
	sed "s,$old,$bse," < $xml > $tmp
	xml=$tmp
    fi
    echo -e "\e[33m$n\e[0m input files for \e[32m${what} stage ${stage}\e[0m"

    copy ${xml} ${top} del

    log_msg "" "Submitting merging job \e[33m${jdl}"
    if test $noact -lt 1 ; then 
	alien_submit alien:${dir}/${jdl} ${run} ${stage} ${tag} ${what}
    fi 
    log_end "" $?
}

# --- Determine the next stage ---------------------------------------
progress()
{
    local out=$1 

    log_msg "" "Deduce next stage for \e[33m$out"
    # First, check for final merge result 
    # echo -e "\nCheck of ${out}/QA_merge_archive.zip"
    alien_ls ${out}/QA_merge_archive.zip > /dev/null 2>&1 
    if test $? -eq 0 ; then 
	stage=6
    else
	#  Then check for production data 
	# echo -e "\nCheck of ${out}/001"
	alien_ls ${out}/001 > /dev/null 2>&1 
	ret=$?
	# echo "ret=$ret"
	if test $ret -ne 0 ; then 
	    # echo "No output, stage 0 to be done"
	    stage=0
	else
	    # Finally, check each merge stage 
	    tmp=0
	    stage=0
	    for i in 4 3 2 1; do 
		# echo -e "\nCheck of ${out}/QA_Stage_${i}"
		alien_ls ${out}/QA_Stage_${i} > /dev/null 2>&1 
		if test $? -ne 0 ; then 
		    tmp=$i 
		else 
		    break
		fi
	    done 
	    stage=$tmp
	fi
    fi 
    log_msg "" "\e[34;m$stage"
    log_end "" 0
}
# --- Upload files ---------------------------------------------------
push()
{
    local bin=$1 
    local data=$2 
    local out=$3
    local tmpdir=`mktemp -d` 

    rm cp.log

    jdls="Run.jdl Merge.jdl Final.jdl"
    for i in $jdls ; do 
	log_msg "" "Creating \e[33m${i}"
	sed -e "s|@out@|${out}|"		\
	    -e "s|@data@|${data}|"		\
	    -e "s|@bin@|${bin}|"		\
	    -e "s|@aliroot@|${aliroot}|"	\
	    -e "s|@root@|${root}|"		\
	    -e "s|@geant@|${geant}|"		\
	    < ${i}.in > ${tmpdir}/${i}
	log_end "" $? 
    done

    log_msg cp.log "Removing and re-creating \e[33m${data}"
    if test $noact -lt 1 ; then 
	alien_rmdir ${data} >> cp.log 2>&1
	alien_mkdir -p $data >> cp.log 2>&1
    fi
    log_end cp.log $?

    local del=1
    if test "X$bin" = "X$data" ; then del=0 ; fi 
    copy run.sh $bin $del
    copy merge.sh $bin $del


    files="simrun.sh		\
	GRP.C			\
	Simulate.C		\
	Config.C		\
	Reconstruct.C		\
	Check.C			\
	Tag.C			\
	QA.C	        	\
	QAConfig.C	       	\
	AOD.C			\
	AODConfig.C		\
	${tmpdir}/Run.jdl	\
	${tmpdir}/Merge.jdl	\
	${tmpdir}/Final.jdl	\
	fmd_corrections.root"
    
    for i in $files ; do 
	copy $i ${data}
    done 
}

# --- Get package versions -------------------------------------------
getVersions()
{
    local ali=$1
    local roo=$2
    local gean=$3

    log_msg "" "Checking software packages"
    if test "x$ali" = "x" ; then 
	log_end "" 1
	log_err "Check versions" "No AliROOT Release specified"
	exit 1
    fi
    if test "x$roo" != "x" && test "x$gean" = "x" ; then 
	log_msg "" "\e[33mAliROOT=$ali ROOT=$roo GEANT=$gean"
	log_end "" 0
	return
    fi
	
    l=`wget -q http://alimonitor.cern.ch/packages/ -O - | \
	sed -n -e '/<tr/,/<\/tr>/ p' | \
	sed -n "/<a.*VO_ALICE@AliRoot::${aliroot}/,/VO_ALICE@ROOT::/ p" | \
	sed -n -e 's/.*VO_ALICE@\(GEAN\|ROO\)T3*::\(v[-0-9a-zA-Z]*\).*/\L\1=\2/gp'|\
	tr '\n' ' '` 
    eval $l
    if test "X$roo" = "X" || test "X$gean" = "X" ; then 
	log_end "" 1 
	log_err "Check versions", "Failed to extract ROOT/GEANT3 versions"
	exit 1
    fi
    root=$roo
    geant=$gean
    log_msg "" "\e[33mAliROOT=$ali ROOT=$root GEANT=$geant"
    log_end "" 0
}

# --- Create an arcive for upload ------------------------------------
archive()
{
    log_msg "" "Creating archive of files" 
    local name=sim_files${version}
    mkdir -p ${name}
    files="\
	run.sh		\
	AOD.C		\
	Check.C		\
	Config.C	\
	doit.sh		\
	Final.jdl.in	\
	GRP.C		\
	Merge.jdl.in	\
	QA.C		\
	README.md	\
	Reconstruct.C	\
	Run.jdl.in	\
	simrun.sh	\
	Simulate.C	\
	Tag.C		\
	merge.sh	\
	fmd_corrections.root"

    for i in $files ; do 
	cp $i ${name}/$i 
    done
    tar -czf ${name}.tar.gz ${name} 
    ret=$?
    rm -rf ${name}
    log_end "" $ret
}

# --- Set some variables ---------------------------------------------
auid=`alien_whoami | sed 's/^ *//'` 
ahome=/alice/cern.ch/user/`echo $auid | sed 's/^\(.\).*/\1/'`/$auid
adir=${ahome}/mc
abin=${ahome}/mc
aout=${ahome}/test

# --- Proces command line options ------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help) usage; exit 0;;
	-t|--tag)       tag=$2    	; shift ;; 
	-i|--id)        id=$2     	; shift ;; 
	-R|--run)    	run=$2    	; shift ;; 
	-e|--generator) gen=$2		; shift ;;
	-c|--copy)      upload=1        ;; 
	-n|--jobs)      jobs=$2   	; shift ;;
	-m|--events)    events=$2 	; shift ;;
	-s|--stage)     stage=$2  	; shift ;; 
	-b|--bin)       abin=$2   	; shift ;; 
	-o|--output)    aout=$2   	; shift ;; 
	-d|--data)      adir=$2   	; shift ;; 
	-a|--aliroot)   aliroot=$2	; shift ;;
	-r|--root)      root=$2   	; shift ;; 
	-g|--geant)     geant=$2  	; shift ;; 
	-f|--final)     minmerge=$2 	; shift ;;
	-A|--archive)   archive   	; exit 0 ;;
	-x|--dry-run)   noact=1		;;
	*) log_err "Unknown option" "$1" ;  exit 1  ;;
    esac
    shift 
done 
abin=$adir



# --- May upload only ------------------------------------------------
if test $upload -gt 0 ; then 
    getVersions "$aliroot" "$root" "$geant"
    push ${abin} ${adir} ${aout}
    if test $stage -lt 0 ; then 
	exit 0
    fi 
fi

# --- Inspect options ------------------------------------------------
case $tag in 
    pp)   id=$tag ; run=118506 ;; 
    PbPb) id=$tag ; run=138190 ;; 
    pPb)  id=$tag ; run=195483 ;;
    Pbp)  id=$tag ; run=196433 ;; 
esac
if test "x$id" = "x" ; then 
    log_err "" "No job identifier given" 
    log_end "" 1
    exit 1
fi
if test "x$run" = "x" ; then 
    log_err "" "No run number given" 
    log_end "" 1
    exit 1
fi
case $stage in 
    0|1|2|3|4|5) : ;; 
    6)  progress $aout/$id/$run ;;
    *)  log_err "Invalid stage" "$stage" ; exit 1 ;;
esac
if test $stage -ge 6 ; then 
    log_msg "" "All done"
    log_end "" 0
    exit 0
fi

# --- Either run the job or merge ------------------------------------
if test $stage -le 0 ; then 
    log_msg "" "Removing old ${aout}/${id}/${run}"
    ret=0
    if test $noact -lt 1 ; then 
	alien_rmdir ${aout}/${id}/${run} > /dev/null 2>&1 
	ret=$?
    fi 
    log_end "" $ret

    log_msg "" "Submitting \e[33mRun.jdl\e[0m for \e[34m$id run\e[0m (\e[34m$jobs\e[0m jobs w/\e[34m$events)"
    if test $noact -lt 1 ; then 
	alien_submit alien:${adir}/Run.jdl ${run} ${jobs} ${events} ${id} ${gen}
	ret=$?
    fi 
    log_end "" $ret
    exit $ret
fi

merge QA ${stage} ${adir} ${aout} ${id} ${run}

#
# EOF
#



