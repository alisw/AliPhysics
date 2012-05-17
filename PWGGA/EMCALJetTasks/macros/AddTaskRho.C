// $Id$

AliAnalysisTaskRho* AddTaskRho(
   const char *outfilename    = "AnalysisResults.root",
   const char *nJets          = "Jets",
   const char *nRho           = "Rho",
   const Double_t minPhi      = 0,
   const Double_t maxPhi      = 2 * TMath::Pi(),
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRho", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRho", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Rho_%s", nJets));
  AliAnalysisTaskRho *rhotask = new AliAnalysisTaskRho(name);
  rhotask->SetJetsName(nJets);
  rhotask->SetRhoName(nRho);
  rhotask->SetJetPhi(minPhi,maxPhi);
  rhotask->SetJetEta(minEta,maxEta);
  rhotask->SetAreaCut(minArea);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput (rhotask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *corho = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(rhotask, 1, corho);

  return rhotask;
}
