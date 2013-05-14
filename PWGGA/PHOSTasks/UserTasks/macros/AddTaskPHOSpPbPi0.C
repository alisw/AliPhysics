AliAnalysisTaskSEPHOSpPbPi0* AddTaskPHOSpPbPi0(UInt_t triggerTag = 0, Bool_t isMCtruth=kFALSE, Bool_t isBadMap=kFALSE, Double_t width = 0.)
{
// Creates a task to analysis pi0 in p-Pb collisions with PHOS
// Author: H. Zhu - 05/14/2013

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSpPbPi0", "No analysis manager to connect to.");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskPHOSpPbPi0", "PHOSpPbPi0 task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }
  if (isMCtruth && type.Contains("ESD")) {
    AliMCEventHandler *mcH = mgr->GetMCtruthEventHandler();
    if (!mcH) {
      ::Error("AddTaskPHOSpPbPi0", "PHOSpPbPi0 task needs the manager to have an MC evnet handler.");
      return NULL;
    }
  }

  const Int_t nBins = 6;
  Float_t cBin[nBins+1] = {0., 20., 40., 60., 80., 90., 100.}; TArrayF tCent(nBins+1, cBin);
  Int_t   buffer[nBins] = { 40,  80,  80,  100, 100, 100    }; TArrayI tBuffer(nBins, buffer);

  Double_t cutsEvent[4] = {  1.,    // 0, min of VtxNCont
                            10.,    // 1, max of VtxZ
                             0.,    // 2, min of Centrality
                           100.  }; // 3, max of Centrality

  Double_t cutsCaloCl[5] = { 0.3,    // 0, min of cluster Energy
                             2.,     // 1, min of NCells
                             0.2,    // 2, min of M02
                             2.5,    // 3, min of DistToBad
                             7e-8 }; // 4, max of TOF

  AliAnalysisTaskSEPHOSpPbPi0 *task = new AliAnalysisTaskSEPHOSpPbPi0("TaskPHOSpPbPi0");
  task->SetUseMC(isMCtruth);
  task->SetXBins(tCent, tBuffer);
  task->SetEventCuts(cutsEvent);
  task->SetCaloClCuts(cutsCaloCl);
  task->SetDecaliWidth(width);    // for decalibration
  task->SetRemovePileup(kFALSE);  // not remove pileup events
  task->SetUseFiducialCut(kTRUE); // use fiducial cut
  task->SetDebugLevel(-1);

  if (triggerTag==0) task->SelectCollisionCandidates(AliVEvent::kINT7);
  if (triggerTag==1) task->SelectCollisionCandidates(AliVEvent::kPHI7);
  if (triggerTag==2) task->SelectCollisionCandidates(AliVEvent::kMB);

  // Bad maps for PHOS
  if (isBadMap) {
    AliOADBContainer badmapContainer("phosBadMap");
    badmapContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
    TObjArray *maps = (TObjArray*)badmapContainer.GetObject(144871,"phosBadMap"); //run number LHC11a

    for (Int_t mod=1;mod<4; mod++) {
      TH2I * h = (TH2I*)maps->At(mod);
      if(h) task->SetPHOSBadMap(mod,h);
    }
  }

  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *output1 = mgr->CreateContainer("histosQA", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  mgr->ConnectOutput(task, 1, output1);
  AliAnalysisDataContainer *output2 = mgr->CreateContainer("histosRD", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  mgr->ConnectOutput(task, 2, output2);

  if (isMCtruth) {
    AliAnalysisDataContainer *output3 = mgr->CreateContainer("histosMC", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
    mgr->ConnectOutput(task, 3, output3);
  }

  return task;
}
