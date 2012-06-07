

AliAnalysisTaskEMCALIsoPhoton *AddTaskEMCALIsoPhoton()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALIsoPhoton", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskEMCALIsoPhoton* ana = new  AliAnalysisTaskEMCALIsoPhoton("");
  
  ana->SelectCollisionCandidates( AliVEvent::kEMC1 | AliVEvent::kMB | AliVEvent::kEMC7 | AliVEvent::kINT7);
  
  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);

  //ana->SetClusThreshold(clusTh);
  ana->SetTrainMode(kTRUE);
  //ana->SetGridMode(kTRUE);
  // ana->SetMcMode(isMC);
  
  AliESDtrackCuts *cutsp = new AliESDtrackCuts;
  cutsp->SetMinNClustersTPC(70);
  cutsp->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  cutsp->SetMaxChi2PerClusterTPC(4);
  cutsp->SetRequireTPCRefit(kTRUE);
  cutsp->SetAcceptKinkDaughters(kFALSE);
  cutsp->SetMaxDCAToVertexZ(3.2);
  cutsp->SetMaxDCAToVertexXY(2.4);
  cutsp->SetDCAToVertex2D(kTRUE);
  cutsp->SetPtRange(0.2);
  cutsp->SetEtaRange(-1.0,1.0);
  ana->SetPrimTrackCuts(cutsp);
  ana->SetPeriod(period.Data());
  if(period.Contains("11"))
    ana->SetGeoName("EMCAL_COMPLETEV1");
  else
    ana->SetGeoName("EMCAL_FIRSTYEARV1");

  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEMCALIsoPhoton", 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
