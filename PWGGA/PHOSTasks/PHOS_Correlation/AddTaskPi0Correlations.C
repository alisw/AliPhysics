AliPHOSCorrelations* AddTaskPi0Correlations (   	const char* name = "Pi0Corr",
							const char* options = "11h",
							const char* fPi0Parametrization = "Wide",
							UInt_t offlineTriggerMask = AliVEvent::kCentral,
							AliPHOSCorrelations::TriggerSelection internalTriggerSelection = AliPHOSCorrelations::kNoSelection,
							Double_t mean = 0.135,
							Double_t sigma = 0.04,
							Int_t downCentLimit = 0,
							Int_t upCentLimit = 90 )
{
	//Author: Ponomarenko Daniil
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
	TString sName = TString (name) + TString (fPi0Parametrization) + centralityBorder ;


	AliPHOSCorrelations* task = new AliPHOSCorrelations(Form("%sTask", sName.Data()),internalTriggerSelection);

	if( TString(options).Contains("10h") )
		task->SetPeriod( AliPHOSCorrelations::kLHC10h );
	if( TString(options).Contains("11h") )
		task->SetPeriod( AliPHOSCorrelations::kLHC11h );
	if( TString(options).Contains("13") )
		task->SetPeriod( AliPHOSCorrelations::kLHC13 );

	// Binning 
	//Any:
	if( AliVEvent::kAny == offlineTriggerMask ) 
	{
		const int nbins = 8;
		Double_t cbin[nbins+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {6, 40, 40, 40, 40, 80, 80, 80};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	//Central:
	if( AliVEvent::kCentral == offlineTriggerMask ) 
	{
		const int nbins = 4;
		Double_t cbin[nbins+1] = {0., 5., 8., 9., 10.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {6, 6, 6, 6};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	// SemiCentral:
	if( AliVEvent::kSemiCentral == offlineTriggerMask ) 
	{
		const int nbins = 8;
		Double_t cbin[nbins+1] = {10., 11., 12., 13., 15., 20., 30., 40., 50.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {40, 40, 40, 40, 40, 40, 40, 40};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	//INT7:
	if( AliVEvent::kINT7 == offlineTriggerMask ) 
	{
		const int nbins = 8;
		Double_t cbin[nbins+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {6, 40, 40, 40, 40, 80, 80, 80};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}
	// MB or PHOS Trigger:
	if( AliVEvent::kMB == offlineTriggerMask || AliVEvent::kPHOSPb == offlineTriggerMask ) 
	{
		const int nbins = 8;
		Double_t cbin[nbins+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
		TArrayD tbin(nbins+1, cbin);
		Int_t    nMixed[nbins] = {6, 40, 40, 40, 40, 80, 80, 80};
		TArrayI tNMixed(nbins, nMixed);
		task->SetCentralityBinning(tbin, tNMixed);
	}


	task->SelectCollisionCandidates(offlineTriggerMask);
	task->SetInternalTriggerSelection(internalTriggerSelection);
	task->EnableTOFCut(true, 100.e-9);
	task->SetMassWindow(mean-sigma, mean+sigma);
	task->SetCentralityBorders(downCentLimit , upCentLimit) ;

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

	TString cname(Form("%sCoutput1", sName.Data()));
	TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), sName.Data()));
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
	mgr->ConnectOutput(task, 1, coutput1);

	return task;
}
