// $Id$

AliHadCorrTask* AddTaskHadCorr(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClusters",
  const char *outClusName    = "CaloClustersCorr",
  const Double_t hadcorr     = 1,
  const Double_t minPt       = 0.15,
  const Double_t phiMatch    = 0.050,
  const Double_t etaMatch    = 0.025,
  const char *outputname     = "AnalysisResults.root",
  const Bool_t   histo       = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHadCorr", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHadCorr", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("HadCorr_%s", outClusName));
  AliHadCorrTask *hcor = new AliHadCorrTask(name, histo);
  hcor->SetTracksName(nTracks);
  hcor->SetClusName(nClusters);
  hcor->SetOutClusName(outClusName);
  hcor->SetPhiMatch(phiMatch);
  hcor->SetEtaMatch(etaMatch);
  hcor->SetHadCorr(hadcorr);
  hcor->SetMinPt(minPt);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(hcor);
    
  // Create containers for input/output
  mgr->ConnectInput (hcor, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cohcor = mgr->CreateContainer(name,
                                                          TList::Class(),
                                                          AliAnalysisManager::kOutputContainer,
                                                          outputname);
  mgr->ConnectOutput(hcor,1,cohcor);
    
  return hcor;
}
