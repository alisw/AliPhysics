detectorQAcontainerName="AD"

runLevelQA()
{
	qaFile=$1
	cp $ALICE_PHYSICS/PWGPP/AD/trending/MakeTrendingADQA.C .
	aliroot -q -b -l .x "MakeTrendingADQA.C(\"${qaFile}\",${runNumber},\"${ocdbStorage}\",kFALSE,kTRUE)"
	#first booleen for Grid connection (true == connection in maketrending) and last boolen for print histo
}

periodLevelQA()
{
	trendingFile=$1
	cp $ALICE_PHYSICS/PWGPP/AD/trending/DrawTrendingADQA.C .
	aliroot -q -b -l .x "DrawTrendingADQA.C(\"${trendingFile}\")"
}
