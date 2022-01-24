#!/bin/bash
### Script has to be executed with source! Otherwise loading the environment does not work as intended

folder=$(pwd)

analysisName=$1
if [[ $1 = "" ]];then
  echo "Please provide name of analysis folder"
fi
analysisFolder=$folder/$analysisName

day=$2
if [[ $2 = "" ]];then
  day=$(date +%Y_%m_%d)
fi
echo "Process systematic from $day"

source $analysisFolder/env.sh

python3 $folder/runCompleteAnalysis.py $analysisFolder $day
