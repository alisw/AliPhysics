AliAnalysisTaskSEPHOSpPbPi0* AddTaskPHOSpPbPi0(Bool_t isMCtruth=kFALSE, UInt_t triggerMask = AliVEvent::kMB, Bool_t isBadMap=kFALSE)
{
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
  Double_t cBin[nBins+1] = {0., 20., 40., 60., 80., 90., 100.}; TArrayD tCent(nBins+1, cBin);
  Int_t    buffer[nBins] = { 40,  80,  80,  100,  100,  100  }; TArrayI tBuffer(nBins, buffer);

  AliAnalysisTaskSEPHOSpPbPi0 *task = new AliAnalysisTaskSEPHOSpPbPi0("TaskPHOSpPbPi0");
  task->SetUseMC(isMCtruth);
  task->SetXBins(tCent, tBuffer);
  task->SetLogWeight(0.06);
  task->SetMinNCells(2);
  task->SetMinClusterEnergy(0.3);
  task->SetMinM02(0.2);
  task->SetMinDistToBad(2.5);
  task->SetDebugLevel(-1);

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

  task->SelectCollisionCandidates(triggerMask); 
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("histosPHOS", TList::Class(), AliAnalysisManager::kOutputContainer, "histosPHOS.root");
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
