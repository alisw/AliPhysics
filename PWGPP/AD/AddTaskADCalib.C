// -*- C++ -*-
// $Id$

AliAnalysisTaskADCalib* AddTaskADCalib(Int_t bcTailMin=16,
				       Int_t bcTailMax=20,
				       Int_t bcExtrapolationMin=10,
				       Int_t bcExtrapolationMax=13) {  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_ADCalib", "No analysis manager found.");
    return MNULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_ADCalib", "This task requires an input event handler");
    return 0;
  }
      
  AliAnalysisTaskADCalib *task = new AliAnalysisTaskADCalib();
  task->SetBCRangeTail(bcTailMin, bcTailMax);
  task->SetBCRangeExtrapolation(bcExtrapolationMin, bcExtrapolationMax);
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("ADCalibListHist",
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s:ADCalib", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput (task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}
