#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber    e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
qaFile=$1
MACRODIR=$ALICE_PHYSICS/PWGPP/EMCAL/QAMacros
SaveImages=1

suffix=".root"
new_suffix=".root"

echo "processing $qaFile" 
echo 

aliroot -b  << EOF
gSystem->AddIncludePath("-I${ALICE_ROOT}/include")
gSystem->SetBuildDir("./",kTRUE)
.L $MACRODIR/CreateEMCALRunQA.C+g
CreateEMCALRunQA("$qaFile","$runNumber","$period","$pass",$SaveImages)
.q
EOF

if [[ $qaFile == *outer* ]] 
then
suffix=_outer.root
new_suffix=_barrel.root
fi

if [[ $qaFile == *barrel* ]]
then
suffix=_barrel.root
new_suffix=_outer.root
fi

if [[ "$suffix" != ".root" ]]
then

mv trending.root trending$suffix
mv ${period}_${pass}_${runNumber}_QAplots.root ${period}_${pass}_${runNumber}_QAplots$suffix

qaFile=${qaFile/$suffix/$new_suffix}

echo
echo "processing $qaFile" 
echo

aliroot -b  << EOF
gSystem->AddIncludePath("-I${ALICE_ROOT}/include")
gSystem->SetBuildDir("./",kTRUE)
.L $MACRODIR/CreateEMCALRunQA.C+g
CreateEMCALRunQA("$qaFile","$runNumber","$period","$pass",$SaveImages)
.q
EOF

mv trending.root trending$new_suffix
mv ${period}_${pass}_${runNumber}_QAplots.root ${period}_${pass}_${runNumber}_QAplots$new_suffix

echo

hadd -v 1 -f trending.root trending_*.root
rm -f trending_*.root

fi

rm -f *.d *.so *.txt

}

periodLevelQA()
{

trendingFile=$1
MACRODIR=$ALICE_PHYSICS/PWGPP/EMCAL/QAMacros
SaveImages=1

echo
echo "Producing PeriodLevel QA plots" 
echo

root -b  << EOF
gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
gSystem->SetBuildDir("./",kTRUE)
.L $MACRODIR/PlotEMCALQATrendingTree.C
PlotEMCALQATrendingTree("$trendingFile",$SaveImages)
.q
EOF

mv trendingPlots.root ${period}_${pass}_trendingPlots.root

echo 

hadd -v 1 -f ${period}_${pass}_EMCALQA.root */*QAplot*.root ${period}_${pass}_trendingPlots.root $trendingFile

echo
echo "Producing QA PDF file" 
echo

root -b -l  << EOF
gSystem->AddIncludePath("-I${ALICE_ROOT}/include")
gSystem->SetBuildDir("./",kTRUE)
.L $MACRODIR/MakeQAPdf.C+g
MakeQAPdf("${period}_${pass}_EMCALQA.root")
.q
EOF

rm -f *.d *.C *.so ${period}_${pass}_trendingPlots.root

}
