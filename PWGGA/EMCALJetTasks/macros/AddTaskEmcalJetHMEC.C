// $Id: AddTaskEmcalJetHMECadron.C 57095 2012-06-20 1:21:07Z mconnors $

AliAnalysisTaskEmcalJetHMEC* AddTaskEmcalJetHMEC(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nJets          = "Jets",
   const char *nTracks        = "PicoTracks",
   const Double_t minPhi      = 1.8,
   const Double_t maxPhi      = 2.74,
   const Double_t minEta      = -0.3,
   const Double_t maxEta      = 0.3,
   const Double_t minArea     = 0.4
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetHMEC", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetHMEC", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Correlations_%s", nJets));
  AliAnalysisTaskEmcalJetHMEC *correlationtask = new AliAnalysisTaskEmcalJetHMEC(name);
  spectratask->SetJetsName(nJets);
  spectratask->SetTracksName(nTracks);
  spectratask->SetJetPhi(minPhi,maxPhi);
  spectratask->SetJetEta(minEta,maxEta);
  spectratask->SetAreaCut(minArea);
 
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTaskMEC(correlationtask);

  // Create containers for input/output
  mgr->ConnectInput (correlationtask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *cojeth = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(correlationtask,1,cojeth);

  return correlationtask;
}
