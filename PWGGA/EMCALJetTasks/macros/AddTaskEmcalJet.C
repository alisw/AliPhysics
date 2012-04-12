AliEmcalJetTask* AddTaskEmcalJet(
						       const char *ntracks        = "Tracks",
						       const char *nclusters      = "CaloClusters",
						       const char *njets          = "Jets",
						       const Int_t a              = 1,
						       const Double_t r           = 0.4,
						       const Int_t t              = 0,
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

  AliEmcalJetTask* jetTask = new AliEmcalJetTask();
  jetTask->SetTracksName(ntracks);
  jetTask->SetClusName(nclusters);
  jetTask->SetJetsName(njets);
  jetTask->SetAlgo(a);
  jetTask->SetMinJetTrackPt(minTrPt);
  jetTask->SetMinJetClusPt(minClPt);
  jetTask->SetRadius(r);
  jetTask->SetType(t);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
  
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 0, coutput1 );
  
  return jetTask;
  
}
