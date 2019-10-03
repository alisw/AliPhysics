detectorQAcontainerName="VZERO"

runLevelQA()
{
	qaFile=$1
	cp $ALICE_PHYSICS/PWGPP/VZERO/trending/MakeTrendingV0QA.C .
	aliroot -q -b -l .x "MakeTrendingV0QA.C(\"${qaFile}\",${runNumber},\"${ocdbStorage}\",kFALSE,kFALSE)"
	#first booleen for Grid connection (true == connection in maketrending) and last boolen for print histo
}

periodLevelQA()
{
	trendingFile=$1
	cp $ALICE_PHYSICS/PWGPP/VZERO/trending/DrawTrendingV0QA.C .
	aliroot -q -b -l .x "DrawTrendingV0QA.C(\"${trendingFile}\")"
}
