#!/bin/bash

nodePrefix=mkrueger_pp_eta_0.80_cutMode
#nodePrefix=mkrueger_pp_mult_100_eta_0.30_cutMode

if [[ $1 == "" ]]; then echo -e "ERROR: No arguments given."; exit; fi

fileDescriptor=$1
inputFolder=$PWD/$1/
if [ ! -e  $inputFolder ]; then echo -e "ERROR: This folder does not exist."; exit; fi

dataFile=${inputFolder}Data.root
mcFile=${inputFolder}MC.root

outputFileCutVariations=${inputFolder}CutVariations.root
if [ -e  $outputFileCutVariations ]; then rm $outputFileCutVariations; fi
outputFileSystematics=${inputFolder}${fileDescriptor}_Sys.root


outputFolderCorr=${inputFolder}Corrected
if [ ! -e  $outputFolderCorr ]; then mkdir $outputFolderCorr; else rm $outputFolderCorr/*; fi

outputFolderPost=${inputFolder}Postprocessed
if [ ! -e  $outputFolderPost ]; then mkdir $outputFolderPost; else rm $outputFolderPost/*; fi

logFileCorr=$outputFolderCorr/log
logFilePost=$outputFolderPost/log
logFileCorrErr=$outputFolderCorr/errlog
logFilePostErr=$outputFolderPost/errlog

for cutSetting in {100..119}
do
	nodeName=${nodePrefix}_$cutSetting
	outputFileCorr=$outputFolderCorr/$cutSetting.root
	aliroot -l -b -q ./applyCorrections.C\(\"$inputFolder\",\"$nodeName\",\"$outputFileCorr\"\) > $logFileCorr 2> $logFileCorrErr

	outputFilePost=$outputFolderPost/$cutSetting.root
	root -l -b -q ../meanPT/meanpt.C\(\"$cutSetting\",\"$outputFolderCorr/\",\"$outputFolderPost/\"\) > $logFilePost 2> $logFilePostErr

	root -l -b -q ./appendHistos.C\(\"$outputFilePost\",\"$cutSetting\",\"$outputFileCutVariations\"\)

done
root -l -b -q ./produceSys.C\(\"$outputFileCutVariations\",\"$outputFileSystematics\"\)
