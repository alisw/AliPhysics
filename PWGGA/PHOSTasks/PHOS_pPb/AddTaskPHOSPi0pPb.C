AliAnalysisTaskPi0Flow* AddTaskPHOSPi0pPb (const char* name = "PHOSPi0pPb",
					   const char* options = "LHC13",
					   const UInt_t offlineTriggerMask = AliVEvent::kINT7,
					   const char* centrality = "V0M",
					   const Int_t nCentBins = 5,
					   const Int_t centEdge0 = 0,
					   const Int_t centEdge1 = 20,
					   const Int_t centEdge2 = 40,
					   const Int_t centEdge3 = 60,
					   const Int_t centEdge4 = 80,
					   const Int_t centEdge5 = 100,
					   const char* trigname = 0)
{
  //Add a task AliAnalysisTaskPi0Flow to the analysis train, for LHC13 PbP data
  //Author: Paul Baetzing
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPi0pPb", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPi0pPb", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPi0Flow* task = new AliAnalysisTaskPi0Flow(Form("%sTask", name));

  if (nCentBins<1) {
    ::Error("AddTaskPHOSPi0pPb", Form("Invalid number of centrality bins: %d",nCentBins));
    return NULL;
  }
  
  Double_t *cbin = new Double_t[nCentBins+1];
  cbin[0] = centEdge0;
  cbin[1] = centEdge1;
  if (nCentBins > 1) cbin[2] = centEdge2;
  if (nCentBins > 2) cbin[3] = centEdge3;
  if (nCentBins > 3) cbin[4] = centEdge4;
  if (nCentBins > 4) cbin[5] = centEdge5;
  TArrayD tbin(nCentBins+1, cbin);
  
  Int_t    *nMixed = new Int_t[nCentBins];
  for (Int_t ibin=0; ibin<nCentBins; ibin++) nMixed[ibin] = 20;
  TArrayI tNMixed(nCentBins, nMixed);
  
  task->SetCentralityBinning(tbin, tNMixed);
  
  if( TString(options).Contains("LHC13") )
    task->SetPeriod( AliAnalysisTaskPi0Flow::kLHC13 );
  
  task->SetCentralityEstimator(centrality);
  
  task->SetEventMixingRPBinning(1);
  task->SelectCollisionCandidates(offlineTriggerMask);
  task->SelectTrigger(trigname);
  task->SetEnablePHOSModule(2, kFALSE);
  task->EnableTOFCut(true, 100.e-9, true);
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
