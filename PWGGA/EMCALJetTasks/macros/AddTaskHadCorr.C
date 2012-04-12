AliHadCorrTask* AddTaskHadCorr(
						       const char *ntracks        = "Tracks",
						       const char *nclusters      = "CaloClusters",
						       const char *outclusname    = "CaloClustersOut",
						       const Double_t hadcorr     = 1,
						       const Double_t minPt       = 0.15
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

  AliHadCorrTask *hcor = new AliHadCorrTask("hcor");
  hcor->SetTracksName(ntracks);
  hcor->SetHadCorr(hadcorr);
  hcor->SetClusName(nclusters);
  hcor->SetMinPt(minPt);
  hcor->SetOutClusName(outclusname);  

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(hcor);
    
  // Create containers for input/output
  mgr->ConnectInput (hcor, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cohcor = mgr->CreateContainer(Form("MatchQAktchemhcorr2"),TList::Class(),AliAnalysisManager::kOutputContainer,Form("Rosi.rho.root"));
  mgr->ConnectOutput(hcor,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(hcor,1,cohcor);
    
  return hcor;
  
}
