AliAnalysisSigmaBarCharged *AddTaskPHOSSigmaBar(
    bool isMC = false, AliVEvent::EOfflineTriggerTypes trig = AliVEvent::kINT7,
    bool isAddHist = false, bool isInvHist = false, bool isQAhist = false,
    Int_t TOFoption = 0, Float_t TOFCut = 150.e-9, Int_t TrackBits = 4,
    Float_t MinECut = 0.6, Int_t NCellsCut = 7, Float_t DispCut = 4.,
    Float_t CPVCut = 10., Float_t TrackEtaCut = 0.8, Int_t TPCClustersCut = 60,
    Float_t TPCsigmasCut = 3., Float_t CPAplusCut = 0.3,
    Float_t CPAminusCut = 0.3, Float_t DCAdaugplusCut = 0.06,
    Float_t DCAdaugminusCut = 0.06, Float_t RADplusCut = 0.25,
    Float_t RADminusCut = 0.15, const char *name = "RR_1") {
  // Add a task AliAnalysisSigmaBarCharged to the analysis train
  // Author: Pavel Gordeev, D.Peresunko, NRC "Kurchatov institute"

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTagging", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTagging", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisSigmaBarCharged *task = new AliAnalysisSigmaBarCharged(name);
  task->SetMC(isMC);
  task->SetAdditionHist(isAddHist);
  task->SetInvMassHist(isInvHist);
  task->SetQAhist(isQAhist);
  if (!isMC) {
    task->SelectCollisionCandidates(trig); // Minimum Bias selection
  }
  task->SetClusterTOF(TOFoption);
  task->SetTOFCut(TOFCut);
  task->SetTrackBits(TrackBits);
  task->SetMinNbarEnergy(MinECut);
  task->SetNcellCut(NCellsCut);
  task->SetDispCut(DispCut);
  task->SetCPVCut(CPVCut);
  task->SetTrackEta(TrackEtaCut);
  task->SetNTPCclusters(TPCClustersCut);
  task->SetTPCsigmas(TPCsigmasCut);
  task->SetCPACut(CPAplusCut, CPAminusCut);
  task->SetDCAdaugCut(DCAdaugplusCut, DCAdaugminusCut);
  task->SetRADCut(RADplusCut, RADminusCut);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      name, TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}
