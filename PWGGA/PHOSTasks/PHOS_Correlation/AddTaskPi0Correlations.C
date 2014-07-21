AliPHOSCorrelations* AddTaskPi0Correlations (   	const char* name = "Pi0Corr",
						const char* options = "11h",
						Double_t sigmaWidth = 3.,
						Int_t downCentLimit = 0,
						Int_t upCentLimit = 90 )
{
	//Author: Ponomarenko Daniil (Daniil.Ponomarenko@cern.ch)
	/* $Id$ */

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) 
	{
		::Error("AddTaskPi0Correlations", "No analysis manager to connect to");
		return NULL;
	}

	if (!mgr->GetInputEventHandler()) 
	{
		::Error("AddTaskPi0Correlations", "This task requires an input event handler");
		return NULL;
	}
	
		
	stringstream ss;
	ss << downCentLimit;
	string strDownCentLimit = ss.str();
	char text[255];
	sprintf(text,"%2i",upCentLimit);
	string strUpCentLimit = text;
	TString centralityBorder = TString ("CB") + strDownCentLimit.c_str() + TString ("t") + strUpCentLimit.c_str() + TString ("Cnt");
	TString sigmaBorder = Form("%2iSigma", int(sigmaWidth*10.));
	if (sigmaWidth == 0) sigmaBorder = "00Sigma";
	TString sName = TString (name) + sigmaBorder + centralityBorder ;


	AliPHOSCorrelations* task = new AliPHOSCorrelations( Form("%sTask", sName.Data()) );

	if( TString(options).Contains("10h") )	{
		task->SetPeriod( AliPHOSCorrelations::kLHC10h );
		task->SetCentralityEstimator("V0M");
	}
	if( TString(options).Contains("11h") )	{
		task->SetPeriod( AliPHOSCorrelations::kLHC11h );
		task->SetCentralityEstimator("V0M");
	}
	if( TString(options).Contains("13") )	{
		task->SetPeriod( AliPHOSCorrelations::kLHC13 );
		task->SetCentralityEstimator("V0A");
	}


	// Binning 
	// TODO: Make other binning for 0-10% and 20-50%
	//Central:
	if( downCentLimit == 0 && upCentLimit == 10 ) 
	{
		const int nbins = 5;
		Double_t cbin[nbins+1] = {0., 2., 4., 6., 8., 10.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {100, 100, 100, 100, 100};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	// SemiCentral:
	if( downCentLimit == 20 && upCentLimit == 50 ) 
	{
		const int nbins = 6;
		Double_t cbin[nbins+1] = {20., 25., 30., 35., 40., 45., 50.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {100, 100, 100, 100, 100, 100};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}

	task->SetAnalysisAlgoritmForReal("ME");
	task->SetAnalysisAlgoritmForMix("ME");
	task->EnableTOFCut(false, 100.e-9);
	task->SelectCollisionCandidates(AliVEvent::kAny);
	task->SetCentralityBorders(downCentLimit , upCentLimit) ;
	task->SetSigmaWidth(sigmaWidth);

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

	TString cname(Form("%sCoutput1", sName.Data()));
	TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), sName.Data()));
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
	mgr->ConnectOutput(task, 1, coutput1);

	return task;
}

