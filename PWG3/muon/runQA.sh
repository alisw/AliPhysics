#!/bin/bash

qaMacroDir="$ALICE_ROOT/PWG3/muon"

baseOutDir=`pwd`


mergeLog="logMerge.txt"
mergeOut=""
terminateDir="${baseOutDir}/terminateRuns"
outTaskName="QAresults.root"

execMerge=1
execMergeAll=1
execTerminate=1
execTrigQA=1
execTrackQA=1
isPrivateProd=0
optList="mateko:pi:"
inputTriggerList=""
while getopts $optList option
do
  case $option in
    m ) execMerge=0;;
    a ) execMergeAll=0;;
    t ) execTerminate=0;;
    e ) execTrigQA=0;;
    k ) execTrackQA=0;;
    o ) outTaskName=$OPTARG;;
    p ) isPrivateProd=1;;
    i ) inputTriggerList=$OPTARG;;
    * ) echo "Unimplemented option chosen."
    EXIT=1
    ;;
  esac
done

shift $(($OPTIND - 1))

# needs 3 arguments
if [[ $# -ne 3 || "$EXIT" -eq 1 ]]; then
    echo "Usage: `basename $0` (-$optList) <runList.txt> <QAx> <alien:///alice/data/20XX/LHCXXy>"
    echo "       -m skip merging (default: run)"
    echo "       -a skip final merging (default: run)"
    echo "       -t skip terminate (default: run)"
    echo "       -e skip run trigger efficiency QA (defult: run)"
    echo "       -k skip run muon QA (defult: run)"
    echo "       -o task output name (default: QAresults.root)"
    echo "       -p is private production. Use directory structure of the plugin"
    echo "       -i input trigger list (default: no list)"
    exit 4
fi

loadLibs="gSystem->Load(\"libANALYSIS.so\");gSystem->Load(\"libOADB.so\");gSystem->Load(\"libANALYSISalice.so\");gSystem->Load(\"libCORRFW.so\");gSystem->Load(\"libPWG3base.so\");"

function mergePerRun()
{
    echo "Merging each run..."
    cd $baseOutDir
    runListName="$1"
    prodDir="$2"
    alienBaseDir="$3"
    aliroot -b <<EOF &> $mergeLog
.L $qaMacroDir/mergeGridFiles.C+
completeProd("${runListName}","${prodDir}","${alienBaseDir}","${outTaskName}",50,"MUON_QA MUON.TriggerEfficiencyMap");
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
${loadLibs}
gSystem->Load("libANALYSIS.so");gSystem->Load("libOADB.so");gSystem->Load("libANALYSISalice.so");gSystem->Load("libCORRFW.so");gSystem->Load("libPWG3base.so");
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
.L $qaMacroDir/terminateQA.C+
CopyFile("$file","$terminateDir",1,"${outTaskName}")
.q 
EOF
	outDir=`grep outDir logCopy.txt | cut -d ":" -f 2`
	#outFile=`grep outFile logCopy.txt | cut -d ":" -f 2`
	forceTerminate=`grep -c "run number not found" logCopy.txt`

	rm logCopy.txt

	cd $outDir
	#ln -s $qaMacroDir/SetAlienHandler.C


	aliroot -b <<EOF &> logTerminate.txt
.x $qaMacroDir/terminateQA.C("${outTaskName}",$forceTerminate)
.q
EOF
	rm logTerminate.txt
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
.x $qaMacroDir/trigEffQA.C+("${runListName}","${outFileName}");
.q
EOF
}

function runTrackQA() {
    echo "Running tracker QA"
    cd $terminateDir
    physSel="$1"
    alienBaseDir="$2"
    lhcPeriod=`echo ${alienBaseDir%"/"} | awk -F "/" ' { print $NF } '`
    aliroot -b <<EOF &> logTrackQA.txt
${loadLibs}
.x $qaMacroDir/PlotMuonQApp.C+("${terminateDir}",0x0,"${inputTriggerList}",${physSel},"${lhcPeriod}","${outTaskName}");
.q
EOF
    cd $baseOutDir
}

# Use absolute path for file inputTriggerList
if [ "${inputTriggerList}" != "" ]; then
  inputTriggerDir=`dirname ${inputTriggerList}`
  if [ "${inputTriggerDir}"="." ]; then
    inputTriggerList="`pwd`/${inputTriggerList}"
  fi
fi

qaProdName="$2"
if [ $isPrivateProd -eq 1 ]; then
    tmpName=${qaProdName//"private"/""}
    if [ "$tmpName" == "$qaProdName" ]; then
	qaProdName="${qaProdName}_private"
    fi
fi

if [ $execMerge -eq 1 ]; then
    mergePerRun $1 $qaProdName $3
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
    minRun=`echo ${mergeOut} | cut -d "_" -f 2`
    maxRun=`echo ${mergeOut} | cut -d "_" -f 3 | cut -d "." -f 1`
    trigOutSuffix=`echo ${qaProdName} | awk -F "/" '{ print $NF }'`
    outName="trigEffQA_${minRun}_${maxRun}_${trigOutSuffix}.root"
    runTrigQA "${mergeOut}" "${outName}"
fi

if [ $execTrackQA -eq 1 ]; then
    runTrackQA 1 $3
fi
