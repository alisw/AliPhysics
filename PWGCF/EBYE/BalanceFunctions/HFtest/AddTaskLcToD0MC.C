//_________________________________________________________//
AliAnalysisTaskLcToD0MC *AddTaskLcToD0MC(Double_t vertexZ=10.,
					 Double_t ptMin=0.0,
					 Double_t ptMax=100.,
					 Double_t etaMin=-0.8,
					 Double_t etaMax=0.8,
					 TString fileNameBase="AnalysisResults") {  
  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLcToD0MC", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskLcToD0MC", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskLcToD0MC *task = new AliAnalysisTaskLcToD0MC("TaskMC");

  // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
  task->SetKinematicsCutsMC(ptMin,ptMax,etaMin,etaMax);

  // vertex cut (x,y,z)
  task->SetVertexDiamond(3.,3.,vertexZ);
  
  mgr->AddTask(task);
    
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer("listQA", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutQA);

  return task;
}
