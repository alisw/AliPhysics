#!/bin/bash

qaMacroDir="$ALICE_PHYSICS/PWGPP/MUON/lite"

baseOutDir=`pwd`


mergeLog="logMerge.txt"
mergeOut=""
terminateDir="${baseOutDir}/terminateRuns"
outTaskName="QAresults.root"
mergeFast=0

execMerge=1
execMergeAll=1
execTerminate=1
execTrigQA=1
execTrackQA=1
execScalers=1
defaultStorage="raw://"
isPrivateProd=0
usePhysicsSelection=1
optList="ad:efi:kmno:pst"
inputTriggerList="0x0"
while getopts $optList option
do
  case $option in
    a ) execMergeAll=0;;
    d ) defaultStorage=$OPTARG;;
    e ) execTrigQA=0;;
    f ) mergeFast=1;;
    i ) inputTriggerList=$OPTARG;;
    k ) execTrackQA=0;;
    m ) execMerge=0;;
    n ) usePhysicsSelection=0;;
    o ) outTaskName=$OPTARG;;
    p ) isPrivateProd=1;;
    s ) execScalers=0;;
    t ) execTerminate=0;;
    * ) echo "Unimplemented option chosen."
    EXIT=1
    ;;
  esac
done

shift $(($OPTIND - 1))

# needs 3 arguments
if [[ $# -ne 3 || "$EXIT" -eq 1 ]]; then
    echo "Usage: `basename $0` (-$optList) <runList.txt> <QAx> <alien:///alice/data/20XX/LHCXXy>"
    echo "       -a skip final merging (default: run)"
    echo "       -d default storage for scalers (default: ${defaultStorage})"
    echo "       -e skip run trigger efficiency QA (default: run)"
    echo "       -f merge fast, skip incomplete prod."
    echo "       -i input trigger list (default: no list)"
    echo "       -k skip run muon QA (defult: run)"
    echo "       -m skip merging (default: run)"
    echo "       -n disable physics selection (default: enable)"
    echo "       -o task output name (default: QAresults.root)"
    echo "       -p is private production. Use directory structure of the plugin"
    echo "       -s skip run scalers trending (default: run)"
    echo "       -t skip terminate (default: run)"
    exit 4
fi

runListName="$1"
qaProdName="$2"
if [ $isPrivateProd -eq 1 ]; then
  tmpName=${qaProdName//"private"/""}
  if [ "$tmpName" == "$qaProdName" ]; then
    qaProdName="${qaProdName}_private"
  fi
fi
alienBaseDir="$3"
lhcPeriod=`echo ${alienBaseDir%"/"} | awk -F "/" ' { print $NF } '`
outFileSuffix="${lhcPeriod}_${qaProdName}"
outFileSuffix=${outFileSuffix//"__"/"_"}

loadAnalysisLibs="gSystem->Load(\"libANALYSIS.so\");gSystem->Load(\"libOADB.so\");gSystem->Load(\"libANALYSISalice.so\");gSystem->Load(\"libCORRFW.so\");gSystem->Load(\"libPWGmuon.so\");gSystem->Load(\"libPWGmuondep.so\");"
includeAliroot="gSystem->AddIncludePath(\"-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include -I${ALICE_INSTALL}/include\");"
includeMuon="gSystem->AddIncludePath(\"-I${ALICE_ROOT}/MUON\");"


function mergePerRun()
{
    echo "Merging each run..."
    cd $baseOutDir

    aliroot -b <<EOF &> $mergeLog
.L $qaMacroDir/mergeGridFiles.C+
completeProd("${runListName}","${qaProdName}","${alienBaseDir}","${outTaskName}",50,"MUON_QA MTR_ChamberEffMap MUON.TrigEfficiencyMap MUON.TriggerEfficiencyMap",${mergeFast});
.q
EOF
}

function mergeRuns()
{
    echo "Merge all runs"
    cd $baseOutDir
    fileListName=$1
    outFilename=$2
    aliroot -b <<EOF &> logMergeAll.txt
${includeAliroot} ${loadAnalysisLibs}
.x $qaMacroDir/mergeGridFiles.C+("${outFilename}","${fileListName}","");
.q
EOF
    
}

function terminateRuns()
{
    echo "Terminating runs..."
    cd $baseOutDir
    currList=`cat $1 | xargs`
    if [ ! -d $terminateDir ]; then
	mkdir $terminateDir
    fi
    for file in $currList; do
	cd $terminateDir
	aliroot -b <<EOF &> logCopy.txt
${includeAliroot} ${loadAnalysisLibs}
.L $qaMacroDir/terminateQA.C+
CopyFile("$file","$terminateDir",1,"${outTaskName}")
.q 
EOF
	outDir=`grep outDir logCopy.txt | cut -d ":" -f 2 | xargs`
	#outFile=`grep outFile logCopy.txt | cut -d ":" -f 2`
	forceTerminate=`grep -c "run number not found" logCopy.txt`

  if [ "`pwd`" != "${outDir}" ]; then
    mv logCopy.txt $outDir/
  fi

	cd $outDir
	#ln -s $qaMacroDir/SetAlienHandler.C


	aliroot -b <<EOF &> logTerminateTrack.txt
${includeAliroot} ${loadAnalysisLibs}
.x $qaMacroDir/terminateQA.C("${outTaskName}",$forceTerminate,0,$usePhysicsSelection)
.q
EOF

aliroot -b <<EOF &> logTerminateTrig.txt
${includeAliroot} ${loadAnalysisLibs}
.x $qaMacroDir/terminateQA.C("${outTaskName}",1,1,$usePhysicsSelection)
.q
EOF
	#rm logTerminate.txt
	#.x $qaMacroDir/runAnalysisTask.C("terminate","grid terminate","",kFALSE,"${outTaskName}");

	#if [ -L "SetAlienHandler.C" ]; then
	#rm SetAlienHandler.C
	#fi
	if [ -e "outputs_valid" ]; then
	    rm outputs_valid
	fi
	cd $baseOutDir
    done
}

function runTrigQA() {
    echo "Running trigger QA"
    cd $baseOutDir
    runListName="$1"
    outFileName="$2"
    aliroot -b <<EOF &> logTrigEffQA.txt
${includeAliroot} ${includeMuon} ${loadAnalysisLibs}
AliLog::SetClassDebugLevel("AliAnalysisTriggerScalers",-1);
.x $qaMacroDir/trigEffQA.C+("${runListName}","${outFileName}","${defaultStorage}",${execScalers});
.q
EOF
}

function runTrackQA() {
    echo "Running tracker QA"
    cd $terminateDir
    physSel="$1"
    aliroot -b <<EOF &> logTrackQA.txt
${includeAliroot} ${loadAnalysisLibs}
.x $qaMacroDir/PlotMuonQA.C+("${terminateDir}",0x0,${inputTriggerList},${physSel},"${outFileSuffix}","${outTaskName}");
.q
EOF
    cd $baseOutDir
}

# Use absolute path for file inputTriggerList
if [ "${inputTriggerList}" != "0x0" ]; then
  inputTriggerDir=`dirname ${inputTriggerList}`
  if [ "${inputTriggerDir}" = "." ]; then
    inputTriggerList="`pwd`/${inputTriggerList}"
  fi
  inputTriggerList="\"${inputTriggerList}\""
fi

if [ $execMerge -eq 1 ]; then
    mergePerRun
fi
mergeOut=`grep -A 1 "Output written" ${mergeLog} | grep -v written`
mergeOutAll=${mergeOut//".txt"/"_merged.txt"}

if [ $execMergeAll -eq 1 ]; then
    mergeRuns $mergeOut "QAresults_Merged.root"
    cp -p $mergeOut $mergeOutAll
    echo "${baseOutDir}/QAresults_Merged.root" >> $mergeOutAll
fi

if [ ! -e $mergeOutAll ]; then
    mergeOutAll=$mergeOut
fi

if [ $execTerminate -eq 1 ]; then
    terminateRuns $mergeOutAll
fi

if [ $execTrigQA -eq 1 ]; then
    outName="trigEffQA_${outFileSuffix}.root"
    runTrigQA "${mergeOut}" "${outName}"
fi

if [ $execTrackQA -eq 1 ]; then
    runTrackQA $usePhysicsSelection
fi
