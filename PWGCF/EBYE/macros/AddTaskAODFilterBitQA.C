//_________________________________________________________//
AliAnalysisTaskAODFilterBitQA *AddTaskAODFilterBitQA(TString taskname = "AODFilterBitQA"
						     ) {
  // Creates an AODFilterBitQA analysis task and adds it to the analysis manager.
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
    ::Error("AddTaskAODFilterBitQA", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (analysisType!="AOD"){
    ::Error("AddTaskAODFilterBitQA", "analysis type NOT AOD --> makes no sense!");
    return NULL;
  }
  
  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskAODFilterBitQA *taskAODFilterBitQA = new AliAnalysisTaskAODFilterBitQA(Form("list_%s",taskname.Data()));

  // ==========================================================================
  // user customization part
  taskAODFilterBitQA->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  // ==========================================================================


  mgr->AddTask(taskAODFilterBitQA);

  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputFilterBitQA";
  AliAnalysisDataContainer *coutAODFilterBitQA = mgr->CreateContainer(Form("list_%s",taskname.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  
  mgr->ConnectInput(taskAODFilterBitQA, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskAODFilterBitQA, 1, coutAODFilterBitQA);

  return taskAODFilterBitQA;
}
