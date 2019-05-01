#!/bin/sh

# Download analysis train results for MC done on pT hard bins and multiple runs
# also merge the runs per pT hard bin into a single file per pT hard bin in the end
# finally extract the histograms into other file and scale each pT hard and merge the final files
#
# how to run
#./DownloadExtractScaleMergePtHardAnalysisFiles.sh RunMinIndex RunMaxIndex pTHardMin pTHardMax path1 path2
#example
#./DownloadExtractScaleMergePtHardAnalysisFiles.sh 0 2 1 21 /alice/sim/2017/LHC17g8a_fast PWGPP/AnalysisQA_AOD/400_20171005-1137
#
# Input files and macros/scripts called:
# - This script expects a list of runs provided in the txt file "runList.txt", 
# each run in different line
# - Also, the merging script "mergePtHardRunFiles.sh" is executed to merge the different 
# histogram files per runs in each pT hard
# - The macro "ScaleExtractPtHardBinHistograms.C" extracts the histograms in 2 separate files
# one scaled with the proper cross section and the other without scaling.
#
#

j=-1

runmin=$1

runmax=$2

pthardmin=$3

pthardmax=$4

path1=$5

path2=$6

# Create pT hard bins directories
for ((bin=$pthardmin;bin < $pthardmax;bin++))
{

mkdir $bin

}

#first loop on run list
for i in `cat runList.txt` 
do

j=`expr $j + 1`

#echo 'Job ' $j  

# Select wich runs from their order in the list
if [   $((j))  -le  $((runmax)) ] 
then

if [   $((j))  -ge  $((runmin)) ] 
then

# Per run loop on the pT hard bin
for ((bin=$pthardmin;bin < $pthardmax;bin++))
{

fullpath=$path1/$bin/$i/$path2/AnalysisResults.root

echo 'pT hard ' $bin  'run: ' $i 'index: ' $j 'copy' $fullpath

alien_cp alien:$fullpath $bin/$i.root

}

fi

fi

done


# Merge the downloaded bins
./mergePtHardRunFiles.sh $pthardmin $pthardmax

# Scale each pT hard bin, extract to file scaled and non scaled histograms
root -q -b -l ScaleExtractPtHardBinHistograms.C


# Merge all the pT hard bins
hadd Scaled.root    */ScaledMerged.root

hadd NotScaled.root */NotScaledMerged.root

























