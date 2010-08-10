#!/bin/sh

# - This script will handle the running of the different Macros of the Mapping of the Conversions for you.
# - In the directory, in which this script is lying in, the different macros and header files which the macros
#   need have to be included as well.
# - Variable #1 is the directory where your rootfiles are put in
# - Variable #2 is the cut you want to analyse, a directory with the same name will be created in #1
# - You need to put in one directory the different root files for MC and data, and rename them so that they 
#   contain the words MC and data an in the corresponding filename. With the outputMerge, outputPdf, outputEps
#   outputSvg and outputGif you can set the different output formats of the files by setting the specific variable
#   to 1.
# - $backup is the the directory where a copy of the macros will be stored.
# - The script will handle all other things for you.

outputMerge=1
outputPdf=0
outputGif=0
outputEps=0
outputSvg=0

cp *.C $1
cp *.h $1

backup= /home/fredi/Dokumente/CERN/ALICE/Photon/Aktuelle\ Makros/newest\ Macros/
cp *.C $backup
cp *.h $backup


cd $1


data=`ls *.root | grep -v Output | grep -i data`
mc=`ls *.root | grep -v Output | grep -i mc`

echo $data
echo $mc

mkdir $2
cp $data $2/
cp $mc $2/
cp *.C $2
cp *.h $2


cd $2

echo $outputMerge 
if [ "$outputMerge" -eq 1 ]; then
	root -l -x -q -b  Photon_Characteristics_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"PhotonCharacteristics_Events_Mult\"\,\"kFALSE\"\,\"gif\"\)
	root -l -x -q -b  Plot_Mapping_Histos_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Mapping_Events_Mult\"\,\"kFALSE\"\,\"gif\"\)
	root -l -x -q -b  Cuts_Events_new.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Cuts_Events_Mult\"\,\"kFALSE\"\,\"gif\"\, \"\"\)
fi

if [ "$outputPdf" -eq 1 ]; then
	mkdir pdf
	root -l -x -q -b Photon_Characteristics_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"PhotonCharacteristics_Events_Mult\"\,\"kTRUE\"\,\"pdf\"\)
	root -l -x -q -b Plot_Mapping_Histos_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Mapping_Events_Mult\"\,\"kTRUE\"\,\"pdf\"\)
	root -l -x -q -b  Cuts_Events_new.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Cuts_Events_Mult\"\,\"kTRUE\"\,\"pdf\"\, \"\"\)
fi

if [ "$outputGif" -eq 1 ]; then
	mkdir gif
	root -l -x -q -b  Photon_Characteristics_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"PhotonCharacteristics_Events_Mult\"\,\"kTRUE\"\,\"gif\"\)
	root -l -x -q -b Plot_Mapping_Histos_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Mapping_Events_Mult\"\,\"kTRUE\"\,\"gif\"\)
	root -l -x -q -b  Cuts_Events_new.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Cuts_Events_Mult\"\,\"kTRUE\"\,\"gif\"\, \"\"\)
fi

if [ "$outputEps" -eq 1 ]; then
	mkdir eps
	root -l -x -q -b  Photon_Characteristics_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"PhotonCharacteristics_Events_Mult\"\,\"kTRUE\"\,\"eps\"\)
	root -l -x -q -b Plot_Mapping_Histos_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Mapping_Events_Mult\"\,\"kTRUE\"\,\"eps\"\)
	root -l -x -q -b  Cuts_Events_new.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Cuts_Events_Mult\"\,\"kTRUE\"\,\"eps\"\, \"\"\)
fi

if [ "$outputSvg" -eq 1 ]; then
	mkdir svg
	root -l -x -q -b  Photon_Characteristics_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"PhotonCharacteristics_Events_Mult\"\,\"kTRUE\"\,\"svg\"\)
	root -l -x -q -b Plot_Mapping_Histos_Events.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Mapping_Events_Mult\"\,\"kTRUE\"\,\"svg\"\)
	root -l -x -q -b  Cuts_Events_new.C\(\"$data\"\,\"$mc\"\,\"$2\"\,\"\"\,\"Cuts_Events_Mult\"\,\"kTRUE\"\,\"svg\"\, \"\"\)
fi

rm $mc
rm $data
rm *.C
rm *.h

cd ..
rm *.C
rm *.h


cd ..
echo Ready with Standard Conversion Analysis

