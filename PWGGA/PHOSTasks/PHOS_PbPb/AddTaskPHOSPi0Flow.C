AliAnalysisTaskPi0Flow* AddTaskPHOSPi0Flow (char* fname="PHOSPi0Flow.root",
					     char* contname="PHOSPi0FlowResults")
{
  //Add a task AliAnalysisTaskPi0Flow to the analysis train
  //Author: Henrik Qvigstad
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPi0Flow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPi0Flow", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPi0Flow* task = new AliAnalysisTaskPi0Flow("PHOSPi0Flow");
  // Reduce binning for reduece memory footprint
  const int kNEdges = 4;
  Double_t cbin[kNEdges] = {0., 10., 40., 80.};
  TArrayD tbin(kNEdges, cbin);
  task->SetCentralityBinning(tbin);
  //task->SetEventMixingRPBinning(9);
  //task->SetMixingArraysLength(10);
  task->SelectCollisionCandidates(AliVEvent::kCentral);
  
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

  // container output into particular file
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,TList::Class(), AliAnalysisManager::kOutputContainer, fname));
  
  return task;
}
