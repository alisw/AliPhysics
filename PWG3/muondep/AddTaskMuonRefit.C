AliAnalysisTaskMuonRefit* AddTaskMuonRefit(Double_t resNB = -1., Double_t resB = -1., Bool_t improveTracks = kFALSE,
					   Double_t sigmaCut = -1., Double_t sigmaCutTrig = -1.)
{
  /// Add AliAnalysisTaskMuonRefit to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonRefit","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskMuonRefit", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonRefit *task = new AliAnalysisTaskMuonRefit("MuonRefit");
  if (!task) {
    Error("AddTaskMuonRefit", "Muon refit task cannot be created!");
    return NULL;
  }
  Double_t clusterResNB[10];
  Double_t clusterResB[10];
  for (Int_t i=0; i<10; i++) { clusterResNB[i] = resNB; clusterResB[i] = resB; }
  task->ResetClusterResolution(clusterResNB, clusterResB);
  task->ImproveTracks(improveTracks, sigmaCut);
  task->SetSigmaCutForTrigger(sigmaCutTrig);
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  return task;
}

