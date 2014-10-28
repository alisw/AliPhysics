

runLevelQA()
{
	qaFile=$1
	detectorQAcontainerName="VZERO"
	cp $ALICE_ROOT/PWGPP/VZERO/trending/MakeTrendingV0QA.C .
	aliroot -q -b -l .x "MakeTrendingV0QA.C(\"${qaFile}\",${runNumber},\"${ocdbStorage}\",kFALSE,KFALSE)"
	#first booleen for Grid connection (true == connection in maketrending) and last boolen for print histo
}

periodLevelQA()
{
	trendingFile=$1
	detectorQAcontainerName="VZERO"
	cp $ALICE_ROOT/PWGPP/VZERO/trending/DrawTrendingV0QA.C .
	aliroot -q -b -l .x "DrawTrendingV0QA.C(\"${trendingFile}\")"
}
