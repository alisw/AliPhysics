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

  // Binning
  // Central:
  if( AliVEvent::kCentral == offlineTriggerMask ) {
    const int nbins = 4;
    Double_t cbin[nbins+1] = {0., 5., 8., 9., 10.};
    TArrayD tbin(nbins+1, cbin);
    Int_t    nMixed[nbins] = {6, 6, 6, 6};
    TArrayI tNMixed(nbins, nMixed);
    task->SetCentralityBinning(tbin, tNMixed);
  }
  // SemiCentral:
  if( AliVEvent::kSemiCentral == offlineTriggerMask ) {
    const int nbins = 8;
    Double_t cbin[nbins+1] = {10., 11., 12., 13., 15., 20., 30., 40., 50.};
    TArrayD tbin(nbins+1, cbin);
    Int_t    nMixed[nbins] = {40, 40, 40, 40, 40, 40, 40, 40};
    TArrayI tNMixed(nbins, nMixed);
    task->SetCentralityBinning(tbin, tNMixed);
  }
  // MB or PHOS Trigger:
  if( AliVEvent::kMB == offlineTriggerMask || AliVEvent::kPHOSPb == offlineTriggerMask ) {
    const int nbins = 8;
    Double_t cbin[nbins+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
    TArrayD tbin(nbins+1, cbin);
    Int_t    nMixed[nbins] = {6, 40, 40, 40, 40, 80, 80, 80};
    TArrayI tNMixed(nbins, nMixed);
    task->SetCentralityBinning(tbin, tNMixed);
  }

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
