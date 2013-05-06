AliAnalysisTaskSEPHOSpPbPi0* AddTaskPHOSpPbPi0(UInt_t triggerTag = 0, Bool_t isMCtruth=kFALSE, Bool_t isBadMap=kFALSE, Double_t logWeight = 0.)
{
// Creates a task to analysis pi0 in p-Pb collisions with PHOS
// H. ZHu - 05/05/2013

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
  Int_t   buffer[nBins] = { 40,  80,  80,  100,  100,  100  }; TArrayI tBuffer(nBins, buffer);

  AliAnalysisTaskSEPHOSpPbPi0 *task = new AliAnalysisTaskSEPHOSpPbPi0("TaskPHOSpPbPi0");
  task->SetUseMC(isMCtruth);
  task->SetXBins(tCent, tBuffer);
  task->SetLogWeight(logWeight); // for ESD decalibration
  task->SetMinDistToBad(2.5);
  task->SetMinNCells(2);
  task->SetMinClusterEnergy(0.3);
  task->SetMinM02(0.2);
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

  AliAnalysisDataContainer *output1 = mgr->CreateContainer("EventInfo", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  AliAnalysisDataContainer *output2 = mgr->CreateContainer("CaloClInfo", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  AliAnalysisDataContainer *output3 = mgr->CreateContainer("Pi0Info", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  mgr->ConnectOutput(task, 1, output1);
  mgr->ConnectOutput(task, 2, output2);
  mgr->ConnectOutput(task, 3, output3);
  if (isMCtruth) {
    AliAnalysisDataContainer *output4 = mgr->CreateContainer("MCInfo", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
    mgr->ConnectOutput(task, 4, output4);
  }

  return task;
}
