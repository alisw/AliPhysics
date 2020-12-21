
AliAnalysisTaskPIDPerformCombIDPtDep *AddTaskPIDCombPerfIDPtDep(TString partType, Int_t fb, Bool_t useCentrality, Float_t centMin, Float_t centMax, TString centEst = "V0M", TString extraString=""){


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  AliAnalysisTaskPIDPerformCombIDPtDep *taskPIDCombMC = new AliAnalysisTaskPIDPerformCombIDPtDep("TaskPIDCombMC");
  taskPIDCombMC->SelectCollisionCandidates(AliVEvent::kINT7);
  taskPIDCombMC->SetFilterBit(fb);

  if (useCentrality){
    taskPIDCombMC->UseCentrality();
    taskPIDCombMC->SetCentralityEstimator(centEst.Data());
    taskPIDCombMC->SetCentralityPercentileRange(centMin, centMax);
    
  }
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("OutputList_%s_%d_%s_%.0f_%.0f_%s", partType.Data(),fb, centEst.Data(), centMin, centMax, extraString.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName.Data());
  
  mgr->AddTask(taskPIDCombMC);
  mgr->ConnectInput(taskPIDCombMC, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPIDCombMC, 1, coutQA);

  return taskPIDCombMC;		     


} 
