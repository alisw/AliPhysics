# 
# Shell script to create a correction/distortion maps
#   
#   tree->histogram
#   histogram->map
#   map->NDfit
#   NDfit->Cheb.regression 
#
#     $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolation.sh
#
makeEnvLocal(){
#
#
# Example usage local at GSI - test farm
#   
    export isBatch=1
    echo makeEnvLocal
    export baliceTPC="ali -n 1"
    ali -n 1
    export isBatch=1
    export batchCommand="qsub -cwd  -V "
}

makeEnvHera(){
#
#
# Example usage GSI big farm
# 
    export isBatch=1
    echo makeEnvHera
    export batchCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=4G  "
    export balice=/hera/alice/$USER/bin/.baliceCalib.sh
    export incScript=$NOTES/aux/rootlogon.C
    source $balice
}




submitTimeDependent(){
    #
    # Input 
    # run.list      - ascii file with runs 
    # residual.list - list of ResidualTree.root
    #  
    # we assume that the run.list and residual list is     
    source $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolation.sh
    source $ALICE_ROOT/libexec/alilog4bash.sh

    # makeEnvLocal; 
    # makeEnvHera; 

    # cp   $ALICE_PHYSICS/../src/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C AliTPCcalibAlignInterpolationMacro.C 
    #
    # 0.) make directory structure
    # 
    #
    alilog_info "BEGIN: 0.) make directory structure"
    wdir=`pwd`
    for arun in `cat  run.list`; do
       mkdir $wdir/$arun;
       cat residual.list | grep $arun > $wdir/$arun/residual.list ;
    done; 
    alilog_info "END: 0.) make directory structure"
    #
    # 1.) Submit query  to get  time dependent info
    #
    alilog_info "BEGIN: 1.) Submit query  to get  time dependent info"
    wdir=`pwd` 
    for arun in `cat run.list`; do
        [ -d $wdir/000$arun ]  &&  	cd $wdir/000$arun;   
        [ -d $wdir/$arun ] && cd $wdir/$arun;       
	echo $wdir `pwd`; 	
        run=`echo $arun| sed s_000__`
	cp $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C AliTPCcalibAlignInterpolationMacro.C
	echo source  $balice >  submitTime.sh
	echo aliroot -b -q $incScript  AliTPCcalibAlignInterpolationMacro.C+\\\(5,$run\\\) >> submitTime.sh
	chmod a+x submitTime.sh  
	if [ $isBatch -eq 0 ] ; then
	    alilog_info "BEGIN:Processing run $arun"
	    ./submitTime.sh  2>&1 | tee  submitTime.log
	    alilog_info "END:Processing run $arun"
        else
	     alilog_info "BEGIN:Processing run $arun on batch"
	    $batchCommand -o submitTime.log submitTime.sh
	     alilog_info "END:Processing run $arun on batch"
	fi;
	cd $wdir
    done    2>&1 | tee  submitScript2.log
    ls $wdir/*/residualInfo.root > timeInfo.list
    alilog_info "END: 1.) Submit query  to get  time dependent info"
    #
    # 2.) Submit hitogramming and filling jobs
    #
    alilog_info "BEGIN: 2.) Submit hitogramming and filling jobs"
    wdir=`pwd` 
    timeDelta=600
    nTracks=100000000;
    for arun in `cat run.list`; do
	wdirRun=$wdir/$arun
	cd $wdirRun
	gminTime=`cat submitTime.log  | grep StatInfo.minTime | gawk '{print $2}' | tail -n 1`	
	gmaxTime=`cat submitTime.log  | grep StatInfo.maxTime | gawk '{print $2}' | tail -n 1`
	rm *.sh
	alilog_info "BEGIN:Processing run $arun on batch"
	for ((itime=$gminTime;itime<$gmaxTime;itime+=$timeDelta));  do   # time loop
	    alilog_info "BEGIN:Processing run $arun at time $itime"
	    echo $itime; 
	    wdirTime=$wdirRun/Time$timeDelta_$itime
            mkdir $wdirTime;
            cd $wdirTime
	    rm -f submit*.sh
	    for ihis in  0 1 2 3 4 5 ; do    #histo time loop
		echo $arun $itime $ihis;
		wdirHis=$wdirTime/his$ihis
		mkdir -p $wdirHis
		cd $wdirHis
		rm -f submit*.sh*
		rm -f *.log
		ln -sf $wdirRun/residual.list .
                ln -sf $wdirRun/residualInfo.root .
		cp $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C AliTPCcalibAlignInterpolationMacro.C
		echo export mapStartTime=$itime > submitHisto$ihis.sh
		echo export mapStopTime=$(($itime+$timeDelta)) >> submitHisto$ihis.sh
		echo source  $balice >>submitHisto$ihis.sh
		echo aliroot -b -q  $incScript AliTPCcalibAlignInterpolationMacro.C+\\\(1,$ihis,$nTracks\\\) >> submitHisto$ihis.sh
		echo cp syswatch.log syswatch_His$ihis.log >>submitHisto$ihis.sh
		echo aliroot -b -q  $incScript AliTPCcalibAlignInterpolationMacro.C+\\\(2,$ihis,0\\\) >> submitHisto$ihis.sh
		echo cp syswatch.log syswatch_Map$ihis.log >>submitHisto$ihis.sh
                chmod a+x submitHisto$ihis.sh
	        echo qsub -cwd  -V -o submitHisto$ihis.log submitHisto$ihis.sh
		if [ $isBatch -eq 0 ] ; then
		    submitHisto$ihis.sh 2>&1 | tee  submitHisto$ihis.log
		else
		    $batchCommand  -o submitHisto$ihis.log submitHisto$ihis.sh
		fi;
		cd $wdirTime	
	    done;
	    alilog_info "END:Processing run $arun at time $itime"
	    cd $wdirRun;
	done;
	alilog_info "END:Processing run $arun"
	cd $wdir;
    done  2>&1 | tee  submitScript2.log
    find $wdir -iname "ResidualMapFull*" > map.list
    #
    # 3.) Submit NDlocalregrression jobs
    #
    alilog_info "BEGIN: 3.) Submit NDlocalregrression jobs"
    wdir=`pwd` 
    find $wdir -iname "ResidualMapFull*"  | grep his1 > map.list
    secStep=2;
    treeNames=( 'deltaRPhiTPCITS' 'deltaRPhiTPCITSTRDDist'  'deltaRPhiTPCITSTOFDist'  'deltaZTPCITSDist'  'deltaZTPCITSTRDDist' 'deltaZTPCITSTOFDist' );
    for amap in `cat map.list`; do
	mapDir=`dirname $amap`; 
	inputFile=`basename $amap`
	ctype=`basename $mapDir | sed s_his__`
	alilog_info "BEGIN:Processing ND fits at directory $mapDir, file $inputFile, type $ctype "
	cd $mapDir;
	for side in  0 1; do
	    for (( sec=0; sec<18; sec+=$secStep )); do 
		theta0=$[side-1];
		theta1=$[side];
		sec0=$sec;
		sec1=$[sec+secStep]
		echo $arun $theta0 $theta1 $sec0 $sec1
		submitScript=submit_${ctype}_${side}_${sec0}
		rm -f $submitScript.*
		echo export inputFile=$inputFile >  ${submitScript}
		echo export inputTree=${treeNames[$ctype]} >> ${submitScript}
		echo export varTheta0=$theta0       >> ${submitScript}
		echo export varTheta1=$theta1       >> ${submitScript}
		echo export varSec0=$sec0           >> ${submitScript}
		echo export varSec1=$sec1           >> ${submitScript}
		echo export runNumber=$arun          >> ${submitScript}
		echo source  $balice >> ${submitScript}
		#cp $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C AliTPCcalibAlignInterpolationMacro.C
		#echo aliroot -b -q  $incScript AliTPCcalibAlignInterpolationMacro.C+\\\(3\\\) >>  ${submitScript}
		echo aliroot -b -q  $incScript $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C+\\\(3\\\) >>  ${submitScript}
		mv ${submitScript} ${submitScript}.sh
		chmod a+x  ${submitScript}.sh    
		alilog_info "BEGIN:Processing ND fits at directory $mapDir, file $inputFile, type $ctype side $side sec = $sec"
		$batchCommand  -o ${submitScript}.log -e ${submitScript}.err  `pwd`/${submitScript}.sh              
	    done;
	done;
	alilog_info "END: Processing ND fits at directory $mapDir, file $inputFile, type $ctype "
	cd $wdir
    done;
    

}


submitProfileHis(){
    #
    # This macros I used to submit callgrind/valgrind jobs
    #    expect inputDir calgrindCommand and calgrindCommand defined
    #  
    # GSI SETUP:
    #
    # inputDir=$NOTES/SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15l.2011test/000240194/Time1445891400/his1
    # calgrindCommand="/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes aliroot "
    # valgrindCommand="/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt -v aliroot"
    #
    wdir=`pwd`
    #
    mkdir $wdir/calgrind/
    cd $wdir/calgrind/
    rm *
    cp $inputDir/*.{list,log,sh,C} .
    cat submitHisto1.sh | sed s_aliroot_"$calgrindCommand"_  | sed s_000_0_ | grep -v "(2,"  >calgrindHisto.sh
    chmod u+x calgrindHisto.sh
    $batchCommand -o calgrindHisto.log -e calgrindHisto.err  calgrindHisto.sh
    #
    mkdir $wdir/valgrind/
    cd $wdir/valgrind/
    rm *
    cp $inputDir/*.{list,log,sh,C} .
    cat submitHisto1.sh | sed s_aliroot_"$valgrindCommand"_  | sed s_000_0_ | grep -v "(2,"> valgrindHisto.sh
    $batchCommand -o valgrindHisto.log -e valgrindHisto.err  valgrindHisto.sh
    mkdir $wdir/calgrind/
    #
    mkdir $wdir/calgrindMap/
    cd $wdir/calgrindMap/
    rm *
    cp $inputDir/*.{list,log,sh,C,root} .
    cat submitHisto1.sh | sed s_aliroot_"$calgrindCommand"_  | sed s_000_0_ | grep -v "(1,"  >calgrindMap.sh
    chmod u+x calgrindMap.sh
    $batchCommand -o calgrindHisto.log -e calgrindHisto.err  calgrindMap.sh
    #
    mkdir $wdir/valgrindMap/
    cd $wdir/valgrindMap/
    rm *
    cp $inputDir/*.{list,log,sh,C,root} .
    cat submitHisto1.sh | sed s_aliroot_"$valgrindCommand"_  | sed s_000_0_ | grep -v "(1,"> valgrindMap.sh
    $batchCommand -o valgrindHisto.log -e valgrindHisto.err  valgrindMap.sh




}

