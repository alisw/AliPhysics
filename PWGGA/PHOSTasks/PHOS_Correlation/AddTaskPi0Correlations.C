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
	
	
	TString className = name;
	TString sigmaName = Form( "%2iSigma", int(sigmaWidth*10.) ) ;
	if( sigmaWidth==0 ) sigmaName = "00Sigma";
	TString sName = Form("%s%sCB%it%iCnt", className.Data(), sigmaName.Data(), downCentLimit, upCentLimit);
	//TString sName = "LOL";

	AliPHOSCorrelations* task = new AliPHOSCorrelations( Form("%sTask", sName.Data()) );

	if( TString(options).Contains("10h") )	
	{
		task->SetPeriod( AliPHOSCorrelations::kLHC10h );
		task->SetCentralityEstimator("V0M");
	}
	
	if( TString(options).Contains("11h") )	
	{
		task->SetPeriod( AliPHOSCorrelations::kLHC11h );
		task->SetCentralityEstimator("V0M");
		if( downCentLimit == 0 && upCentLimit == 10 ) 
		{
			Double_t meanParametrs[2]  = {-0.000129767, 0.138874 };
			Double_t sigmaParametrs[4] = {5.73226e-06, -0.00879368, 0.00462739 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}

		if( downCentLimit == 20 && upCentLimit == 50 ) 
		{
			Double_t meanParametrs[2]  = {-8.35555e-05, 0.136538 };
			Double_t sigmaParametrs[4] = {-7.61949e-06, 1.20701e-06, 0.00474992 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
	}
	
	if( TString(options).Contains("13") )	
	{
		task->SetPeriod( AliPHOSCorrelations::kLHC13 );
		task->SetCentralityEstimator("V0A");
		if( downCentLimit == 0 && upCentLimit == 10 ) 
		{
			Double_t meanParametrs[2]  = {-4.64539e-05, 0.134773 };
			Double_t sigmaParametrs[3] = {0.00383029, 0.0041709, 0.00468736 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}

		if( downCentLimit == 20 && upCentLimit == 50 ) 
		{
			Double_t meanParametrs[2]  = {-4.90799e-06, 0.134566 };
			Double_t sigmaParametrs[4] = {0.00293721, 0.00622308, 0.00468625 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
	}


	// Mixed binning 
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

	// Events
	task->SelectCollisionCandidates(AliVEvent::kAny);
	task->SetCentralityBorders((Double_t)downCentLimit , (Double_t)upCentLimit) ;
	// Clasters
	task->EnableTOFCut(false, 100.e-9);
	task->SwitchOnPionEfficiency();
	task->SwitchOnMassParametrisation();
	task->SetSigmaWidth(sigmaWidth);
	// Tracks
	task->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    task->SetTrackStatus(AliVTrack::kITSrefit);
    task->SetTPCSharedClusterFraction(0.4);
    task->SwitchOnAODTrackSharedClusterSelection();
    task->SwitchOffTrackHitSPDSelection();
    task->SetTrackFilterMask(786);


	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

	TString cname(Form("%sCoutput1", sName.Data()));
	TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), sName.Data()));
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
	mgr->ConnectOutput(task, 1, coutput1);

	return task;
}

