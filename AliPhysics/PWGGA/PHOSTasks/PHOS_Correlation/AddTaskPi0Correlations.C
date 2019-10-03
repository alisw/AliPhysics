AliPHOSCorrelations* AddTaskPi0Correlations (   	const char* name = "Pi0Corr",
						const char* options = "11h",
						Int_t downCentLimit = 0,
						Int_t upCentLimit = 90,
						const char* suffix = "" )
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
	TString sName = Form("%s%sCB%it%iCnt", className.Data(), suffix, downCentLimit, upCentLimit);

	TString combinedName;
	combinedName.Form("%sTask", sName.Data());

	AliPHOSCorrelations* task = new AliPHOSCorrelations( combinedName );

	// Mixing binning 
	if( downCentLimit == 0 && upCentLimit == 10 ) 
	{
		//Central
		const int nbins = 5;
		Double_t cbin[nbins+1] = {0., 2., 4., 6., 8., 10.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
		task->SetTriggerSelectionInPbPb(AliPHOSCorrelations::kCent);
	}
	else
	if( downCentLimit == 0 && upCentLimit == 20 ) 
	{
		//Central
		const int nbins = 4;
		Double_t cbin[nbins+1] = {0., 5., 10., 15., 20.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
		task->SetTriggerSelectionInPbPb(AliPHOSCorrelations::kCentAndSCent);
	}
	else
	if( downCentLimit == 10 && upCentLimit == 20 ) 
	{
		//Central
		const int nbins = 5;
		Double_t cbin[nbins+1] = {10., 12., 14., 16., 18., 20.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
		task->SetTriggerSelectionInPbPb(AliPHOSCorrelations::kSCent);
	}
	else
	if( downCentLimit == 20 && upCentLimit == 40 ) 
	{
		// SemiCentral
		const int nbins = 4;
		Double_t cbin[nbins+1] = {20., 25., 30., 35., 40.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10,};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	else
	if( downCentLimit == 20 && upCentLimit == 50 ) 
	{
		// SemiCentral
		const int nbins = 6;
		Double_t cbin[nbins+1] = {20., 25., 30., 35., 40., 45., 50.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	else
	if( downCentLimit == 40 && upCentLimit == 60 ) 
	{
		// SemiCentral
		const int nbins = 4;
		Double_t cbin[nbins+1] = {40., 45., 50., 55., 60.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	else
	if( downCentLimit == 60 && upCentLimit == 90 ) 
	{
		// SemiCentral
		const int nbins = 6;
		Double_t cbin[nbins+1] = {60., 65., 70., 75., 80., 85., 90.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {10, 10, 10, 10, 10, 10};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	else
	{
		::Error("AddTaskPi0Correlations", Form("No centrality setings for %i-%i\%. Skip", downCentLimit, upCentLimit) );
		return NULL;
	}


	// Period depended properties
	if( TString(options).Contains("11h") )	
	{
		task->SetCentralityEstimator("V0M");

		if( downCentLimit == 0 && upCentLimit == 10 ) 
		{
			Double_t meanParametrs[2]  = {-3.36111e-05, 0.137305 };
			Double_t sigmaParametrs[3] = {0.0081226, -3.63246e-08, 0.00229115 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
		else
		{
			::Error("AddTaskPi0Correlations", Form("No mass window setings for %i-%i\% in %s. Set default", downCentLimit, upCentLimit, options));
			Double_t meanParametrs[2]  = {-2.72612e-05, 0.1357 };
			Double_t sigmaParametrs[3] = {0.00499517, 0.00950036, 0.00275293 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
	}
	else
	if( TString(options).Contains("13") )	
	{
		task->SetCentralityEstimator("V0A");

		if( downCentLimit == 0 && upCentLimit == 10 ) 
		{
			Double_t meanParametrs[2]  = {-1.14244e-05, 0.134433 };
			Double_t sigmaParametrs[3] = {-0.00525797, 0.00212986, 0.00407932 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
		else
		{
			::Error("AddTaskPi0Correlations", Form("No mass window setings for %i-%i\% in %s. Set default", downCentLimit, upCentLimit, options));
			Double_t meanParametrs[2]  = {2.40957e-05, 0.13428 };
			Double_t sigmaParametrs[3] = {-0.00293592, 0.00624031, 0.0042106 };
			task->SetMassMeanParametrs(meanParametrs);
			task->SetMassSigmaParametrs(sigmaParametrs);
		}
	}
	else
	{
		::Error("AddTaskPi0Correlations", Form("Unknown period \"%s\". Skip", options));
		return NULL;
	}

	// Period setup
	task->SetPeriodName( TString(options) );
	// Events
	task->SelectCollisionCandidates(AliVEvent::kAny);
	task->SetCentralityBorders((Double_t)downCentLimit , (Double_t)upCentLimit) ;
	task->SetEventPlaneMethod("V0");
	// Mixing
	task->SetEventMixingRPBinning(5);
	task->SetEventMixingVtxBinning(5);
	// Clasters
	task->EnableTOFCut(false, 100.e-9);
	task->SwitchOnPionEfficiency();
	task->SwitchOnMassParametrisation();
	task->SetSigmaWidth(3.);
	// Tracks
	task->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    task->SetTrackStatus(AliVTrack::kITSrefit);
    task->SetTPCSharedClusterFraction(0.4);
    task->SwitchOnAODTrackSharedClusterSelection();
    task->SwitchOffTrackHitSPDSelection();
    task->SetTrackFilterMask(786);

    task->ShowTaskInfo();
	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

	TString cname(Form("%sCoutput1", combinedName.Data()));
	TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), combinedName.Data()));
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
	mgr->ConnectOutput(task, 1, coutput1);

	return task;
}

