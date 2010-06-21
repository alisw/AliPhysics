#! /bin/bash
#
#
#
# This script gests as input a directory where the AnalysisResults root file is stored
# , it also needs the desired output directory where the produced root files are put.
# If nothing is given it will use ./ for the input directory and ./Output for the output
#
#Input 1: Root file to analyze (not including the .root) Default: AnalyisResults
#Input 2: Input directory  Default:$PWD 
#Input 3: Output directory Default: $PWD/Output  (directory will be created if it does not exist)
#

RootFile="";
if [ -n $1 ]; then
    RootFile=AnalysisResults
else
    RootFile=$1
fi


InputDirectory="";
if [ -n $2 ]; then
    InputDirectory=$PWD/
else
    InputDirectory=$2
fi
echo Input directory is $InputDirectory

OutputDirectory="";
if [ -n $3 ]; then
    OutputDirectory=$PWD/Output/
else
    OutputDirectory=$3
fi
 
if [ ! -d $OutputDirectory ]; then
    mkdir $OutputDirectory
fi 
echo Output directory is $OutputDirectory

DatFilenameBase=RB-data-AnalysisResults;
Suffix=gif;

if [ -f $OutputDirectory/$DatFilenameBase.dat ]; then
    echo Warning: The file RB-data-AnalysisResults.dat exists, please remove before continuing.
    echo -e "\t Otherwise the file will accumulate more entries than it is supposed to."
    exit;
fi


root -b -q $ALICE_ROOT/PWG4/GammaConv/macros/MakeCutLog.C\(\"$RootFile\"\,\"$InputDirectory\"\,\"$OutputDirectory\"\)

# Read the different cuts form the Cut selection log file
exec<"$OutputDirectory/CutSelection.log"

while read cutSelection
do
    echo CutSelection is $cutSelection;
    root -b -q $ALICE_ROOT/PWG4/GammaConv/macros/Extract_IntegratedPi0Yield.C\(\"$cutSelection\"\,\"$RootFile\"\,\"$InputDirectory\"\,\"$OutputDirectory\"\"\);
    
    root -b -q $ALICE_ROOT/PWG4/GammaConv/macros/Extract_Pi0_Characteristics.C\(\"$cutSelection\"\,\"$RootFile\"\,\"$InputDirectory\"\,\"$OutputDirectory\"\);

done

    root -b -q $ALICE_ROOT/PWG4/GammaConv/macros/Plot_IntegratedPi0Yield.C\(\"$DatFilenameBase\"\,\"$OutputDirectory\"\,\"$Suffix\"\)

