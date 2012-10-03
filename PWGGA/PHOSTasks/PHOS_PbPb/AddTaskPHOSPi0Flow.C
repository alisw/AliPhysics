AliAnalysisTaskPi0Flow* AddTaskPHOSPi0Flow (const char* name = "PHOSPi0Flow",
					    const char* options = "",
					    UInt_t offlineTriggerMask = AliVEvent::kCentral )
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

  AliAnalysisTaskPi0Flow* task = new AliAnalysisTaskPi0Flow(Form("%sTask", name));

  // Reduce binning for reduece memory footprint
  const int nbins = 3;
  Double_t cbin[nbins+1] = {0., 10., 40., 80.};
  TArrayD tbin(nbins+1, cbin);
  Int_t    nMixed[nbins] = {4, 20, 50};
  TArrayI tNMixed(nbins, nMixed);
  task->SetCentralityBinning(tbin, tNMixed);

  //task->SetEventMixingRPBinning(9);
  //task->SetMixingArraysLength(10);
  task->SelectCollisionCandidates(offlineTriggerMask);
  
  if( TString(options).Contains("11h") )
    task->SetPeriod( AliAnalysisTaskPi0Flow::kLHC11h );

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
