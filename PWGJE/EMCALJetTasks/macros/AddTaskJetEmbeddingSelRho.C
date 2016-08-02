AliJetEmbeddingSelRhoTask* AddTaskJetEmbeddingSelRho(
	const char *rhoname,
	Double_t    minrho,
	Double_t    maxrho,
	const char     *idname       = "",
	const char     *tracksName   = "Tracks",
	const char     *clusName     = "CaloClustersCorr",
	const Double_t  minPt        = 10,
	const Double_t  maxPt        = 100,
	const Double_t  minEta       = -0.9,
	const Double_t  maxEta       = 0.9,
	const Double_t  minPhi       = 0,
	const Double_t  maxPhi       = TMath::Pi() * 2,
	const Int_t     nTracks      = 1,
	const Int_t     nClus        = 0,
	const Bool_t    copyArray    = kFALSE,
	const char     *subwagonname = ""
	) 
{
	
	
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AliJetEmbeddingSelRhoTask", "No analysis manager to connect to.");
		return NULL;
	}  
	
	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler())
	{
		::Error("AliJetEmbeddingSelRhoTask", "This task requires an input event handler");
		return NULL;
	}
	
	//-------------------------------------------------------
	// Init the task and do settings
	//-------------------------------------------------------
	idname+=subwagonname;
	
	TString taskName = Form("TaskEmb%sSel%s_%.0f_%.0f", idname, rhoname, minrho*10., maxrho*10.);
	
	AliJetEmbeddingSelRhoTask *jetEmb = new AliJetEmbeddingSelRhoTask(taskName);
	
	jetEmb->SetRhoName(rhoname);
	jetEmb->SetRhoRange(minrho, maxrho);
	
	jetEmb->SetTracksName(tracksName);
	jetEmb->SetClusName(clusName);
	jetEmb->SetEtaRange(minEta, maxEta);
	jetEmb->SetPhiRange(minPhi, maxPhi);
	jetEmb->SetPtRange(minPt, maxPt);
	jetEmb->SetCopyArray(copyArray);
	jetEmb->SetNClusters(nClus);
	jetEmb->SetNTracks(nTracks);
	
	//-------------------------------------------------------
	// Final settings, pass to manager and set the containers
	//-------------------------------------------------------
	
	mgr->AddTask(jetEmb);
    
	// Create containers for input/output
	mgr->ConnectInput (jetEmb, 0, mgr->GetCommonInputContainer() );
	
	//Connect output
	TString contName(taskName);
	//contName+="Input";
	TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
	mgr->ConnectOutput(jetEmb,1,coutput1);
	
	contName+="Input";
	TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
	mgr->ConnectOutput(jetEmb,2,coutput2);
	
	return jetEmb;
}

