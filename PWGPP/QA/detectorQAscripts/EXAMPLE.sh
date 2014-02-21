# this is a simple template
# it defines two functions: runLevelQA and periodLevelQA
# 
#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumer     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
  #full path of QAresults.root is provided
  qaFile=$1

  #aliroot....

  #should produce a file trending.root
  touch trending.root
}

periodLevelQA()
{
  #path of the merged period trending.root is provided
  trendingFile=$1 
  #runs in the production dir: ...../LHCXXx/passX/
  #the running dir contains all the run sub directories
  #named like 000123123 with the outputs of runLevelQA

  #aliroot...
}
