#
# shel script to genarate certificate
# arguments:
# 
amacro=$1
esdList=$2
echo Get certificate for $amacro
aliroot -l -b -q  runPWG1Train.C\(\"$amacro\", \"$esdList\"\)  2>&1 | tee out.log
makeSummary.sh > summary.log






