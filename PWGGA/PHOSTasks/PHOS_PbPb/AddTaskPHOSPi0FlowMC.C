AliAnalysisTaskPi0FlowMC* AddTaskPHOSPi0FlowMC (const char* name = "PHOSPi0FlowMC",
					    const char* options = ""
						// , UInt_t offlineTriggerMask = AliVEvent::kCentral
						)
{
  //Add a task AliAnalysisTaskPi0FlowMC to the analysis train
  //Author: Henrik Qvigstad
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
  // const int nbins = 8;
  // Double_t cbin[nbins+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
  // TArrayD tbin(nbins+1, cbin);
  // Int_t    nMixed[nbins] = {6, 40, 40, 40, 40, 80, 80, 80};
  // TArrayI tNMixed(nbins, nMixed);
  // task->SetCentralityBinning(tbin, tNMixed);

  // Binning
  const int nbins = 6;
  Double_t cbin[nbins+1] = {0., 5., 10., 20., 40., 60., 80.};
  TArrayD tbin(nbins+1, cbin);
  Int_t    nMixed[nbins] = {6, 6, 40, 40, 40, 80};
  TArrayI tNMixed(nbins, nMixed);
  task->SetCentralityBinning(tbin, tNMixed);

  //task->SetEventMixingRPBinning(9);
  //task->SetMixingArraysLength(10);
  //task->SelectCollisionCandidates(offlineTriggerMask);
  

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
