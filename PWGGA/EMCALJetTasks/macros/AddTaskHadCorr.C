// $Id$

AliHadCorrTask* AddTaskHadCorr(
  const char *ntracks        = "Tracks",
  const char *nclusters      = "CaloClusters",
  const char *outclusname    = "CaloClustersCorr",
  const Double_t hadcorr     = 1,
  const Double_t minPt       = 0.15,
  const char *outputname     = "HadCorrOutput.root"
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

  TString name(Form("HadronicCorrection_%s", outclusname));
  AliHadCorrTask *hcor = new AliHadCorrTask(name);
  hcor->SetTracksName(ntracks);
  hcor->SetClusName(nclusters);
  hcor->SetOutClusName(outclusname);  
  hcor->SetHadCorr(hadcorr);
  hcor->SetMinPt(minPt);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(hcor);
    
  // Create containers for input/output
  mgr->ConnectInput (hcor, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cohcor = mgr->CreateContainer(name,
                                                          TList::Class(),AliAnalysisManager::kOutputContainer,
                                                          outputname);
  mgr->ConnectOutput(hcor,1,cohcor);
    
  return hcor;
}
