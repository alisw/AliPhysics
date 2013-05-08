AliAnalysisTaskPi0FlowMC* AddTaskPHOSPi0pPbMC (const char* name = "PHOSPi0pPbMC",
					    const char* options = "",
					       UInt_t offlineTriggerMask = AliVEvent::kINT7)
{
  //Add a task AliAnalysisTaskPi0FlowMC to the analysis train, for LHC13 PbP data
  //Author:  Paul Baetzing
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPi0FlowMC", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPi0FlowMC", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPi0FlowMC* task = new AliAnalysisTaskPi0FlowMC(Form("%sTask", name));


  // Binning
  const int nbins = 5;
  Double_t cbin[nbins+1] = {0., 20., 40., 60., 80.,100.};
  TArrayD tbin(nbins+1, cbin);
  Int_t    nMixed[nbins] = {40, 40, 40, 40, 40};
  TArrayI tNMixed(nbins, nMixed);
  task->SetCentralityBinning(tbin, tNMixed);

  //task->SetEventMixingRPBinning(9);
  //task->SetMixingArraysLength(10);
  task->SelectCollisionCandidates(offlineTriggerMask);
  

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
