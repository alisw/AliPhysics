AliAnalysisTaskPi0Flow* AddTaskPHOSPi0Flow ()
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
  const int nbins = 3;
  Double_t cbin[nbins+1] = {0., 10., 40., 80.};
  TArrayD tbin(nbins+1, cbin);
  task->SetCentralityBinning(tbin);
  Int_t    nMixed[nbins] = {4, 20, 50};
  TArrayI tNMixed(nbins, nMixed);
  task->SetNMixedPerCentrality(tNMixed);

  //task->SetEventMixingRPBinning(9);
  //task->SetMixingArraysLength(10);
  task->SelectCollisionCandidates(AliVEvent::kCentral);
  

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PHOSPi0FlowCoutput1", TList::Class(), AliAnalysisManager::kOutputContainer, 
							    Form("%s:PHOSPi0Flow", AliAnalysisManager::GetCommonFileName()) 		);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
