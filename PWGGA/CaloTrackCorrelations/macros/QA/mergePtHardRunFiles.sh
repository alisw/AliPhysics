
#!/bin/sh

# how to run
#./mergePtHardRunFiles.sh  pTHardMin pTHardMax 
#
# example
#./mergePtHardRunFiles.sh  1 21 
#
# This script is called by  DownloadExtractScaleMergePtHardAnalysisFiles.sh
#

pthardmin=$1

pthardmax=$2

for ((bin=$pthardmin;bin < $pthardmax;bin++))
{

echo 'Merge pT hard ' $bin  

hadd $bin/Merged.root $bin/*root

}

