//_________________________________________________________//
AliAnalysisTaskParticleStudies *AddTaskParticleStudies(TString taskname = "ParticleStudies", 
						       Double_t etaMin = -0.8,
						       Double_t etaMax = 0.8,
						       Double_t ptMin = 0.13,
						       Double_t ptMax = 4.0
						       ) {
  // Creates an ParticleStudies analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBalancePsiCentralityTrain", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskParticleStudies", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskParticleStudies *taskParticleStudies = new AliAnalysisTaskParticleStudies(Form("list_%s",taskname.Data()));

   // ==========================================================================
  // user customization part
  taskParticleStudies->SetEtaRange(etaMin,etaMax);
  taskParticleStudies->SetPtRange(ptMin,ptMax);
  // ==========================================================================

  mgr->AddTask(taskParticleStudies);

  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.ParticleStudies";
  AliAnalysisDataContainer *coutParticleStudies = mgr->CreateContainer(Form("list_%s",taskname.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  
  mgr->ConnectInput(taskParticleStudies, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskParticleStudies, 1, coutParticleStudies);

  return taskParticleStudies;
}
