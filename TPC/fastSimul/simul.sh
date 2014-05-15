#!/bin/sh
main()
{
  #
  # run in proper action depending on the selection
  #  
  if [[ $# -lt 1 ]]; then
    if [[ ! "$0" =~ "bash" ]]; then
      echo " Please select action"
    fi
    return
  fi
  runMode=$1
  umask 0002
  shift
  case $runMode in
   "runJob") runJob "$@";;
    "makeEnvLocal") makeEnvLocal "$@";;
    "makeSubmitRun") makeSubmitRun "$@";;
   *) 
   eval "${runMode} $@" 
   ;;
  esac
  return;
}



exampleCase(){
#
#  Example case to subit Toy MC jobs
# 
   source $ALICE_ROOT/TPC/fastSimul/simul.sh
   makeEnvLocal
   makeSubmitRUN 80 400
   ls `pwd`/MC*/trackerSimul.root >  trackerSimul.list

}



runJob()
{
#runFastMCJob      
    echo  $ROOTSYS
    which root.exe
    which aliroot
    echo PWD `pwd`
    ntracks=$1
    echo Submitting ntracks = $ntracks
    echo command aliroot  -q -b  "$mcPath/simul.C\($ntracks\)"    
    command aliroot  -q -b  "$mcPath/simul.C($ntracks)"    
    return;
}


makeEnvLocal(){
#
#
# Example usage local 
# jobs to be submitted form the lxb1001 or lxb1002
#(here we have 80 nodes and user disk)
# 
    echo makeEnvLocal
    export baliceTPC=$HOME/.baliceTPC
    export mcPath=$ALICE_ROOT/TPC/fastSimul
    export batchCommand="qsub -cwd  -V "
}

makeSubmitRUN(){
#
# submits jobs
#   
    wdir=`pwd`;
    njobs=$1
    ntracks=$2
    for (( job=1; job <= $njobs; job++ ));  do  
	echo $job;  
	mkdir $wdir/MC$job
	cd $wdir/MC$job
 	echo $batchCommand    -o  toyMC.log  $mcPath/simul.sh runJob  $ntracks
 	$batchCommand    -o  toyMC.log  $mcPath/simul.sh runJob $ntracks
	cd $wdir
    done 
}



main "$@"

