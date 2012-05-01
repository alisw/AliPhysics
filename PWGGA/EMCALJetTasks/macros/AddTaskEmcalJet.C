// $Id$ 

AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClusters",
  const char *nJets          = "Jets",
  const Int_t algo           = 1,
  const Double_t radius      = 0.4,
  const Int_t type           = 0,
  const Double_t minTrPt     = 0.15,
  const Double_t minClPt     = 0.15
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskAliEmcalJet", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskAliEmcalJet", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("JetTask_%s", nJets));
  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  jetTask->SetTracksName(nTracks);
  jetTask->SetClusName(nClusters);
  jetTask->SetJetsName(nJets);
  jetTask->SetAlgo(algo);
  jetTask->SetMinJetTrackPt(minTrPt);
  jetTask->SetMinJetClusPt(minClPt);
  jetTask->SetRadius(radius);
  jetTask->SetType(type);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (jetTask, 0, cinput);

  return jetTask;
}
