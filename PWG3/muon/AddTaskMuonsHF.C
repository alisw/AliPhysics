AliAnalysisTaskSEMuonsHF* AddTaskMuonsHF(Int_t mode=0, Bool_t isMC=kFALSE, Bool_t isTree=kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuonsHF", "No analysis manager to connect to.");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskMuonsHF", "MuonsHF task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }
  if (isMC && type.Contains("ESD")) {
    AliMCEventHandler *mcH = mgr->GetMCtruthEventHandler();
    if (!mcH) {
      ::Error("AddTaskMuonsHF", "MuonsHF task needs the manager to have an MC evnet handler.");
      return NULL;
    }
  }
  if (isTree) {
    AliAODHandler *aodH = (AliAODHandler*)mgr->GetOutputEventHandler();
    if (!aodH) {
      ::Error("AddTaskMuonsHF", "MuonsHF task needs the manager to have an AOD output handler.");
      return NULL;
    }
  }

  // set cuts for events or muons selection
  Double_t cutsEvsH[3] ={      0.5,   // low limit of Ncontrs
                          999999.0,   // up limit of |vz|
                          999999.0};  // up limit of vt
  Double_t cutsMuon[12]={-999999.0,   //   0, min of 3-momentum
                          999999.0,   //   1, max of 3-momnentum
                         -999999.0,   //   2, PtMin
                          999999.0,   //   3, PtMax
                              -4.0,   //   4, EtaMin
                              -2.5,   //   5, EtaMax
                         -999999.0,   //   6, DCAmin
                          999999.0,   //   7, DCAmax
                               1.5,   //   8, for trigger
                               3.5,   //   9, for trigger
                             171.0,   //  10, ThetaAbsEndMin
                             178.0 };  // 11, ThetaAbsEndMax
  AliAnalysisTaskSEMuonsHF *taskMuonsHF = new AliAnalysisTaskSEMuonsHF("MuonsHF Analysis Task");
  taskMuonsHF->SetAnaMode(mode);
  taskMuonsHF->SetUseMC(isMC);
  taskMuonsHF->SetIsOutputTree(isTree);
  taskMuonsHF->SetEvsHCuts(cutsEvsH);
  taskMuonsHF->SetMuonCuts(cutsMuon);
  taskMuonsHF->SetDimuCuts(cutsMuon);
  taskMuonsHF->SelectCollisionCandidates();
  mgr->AddTask(taskMuonsHF);

  mgr->ConnectInput(taskMuonsHF, 0, mgr->GetCommonInputContainer());
  if (isTree) mgr->ConnectOutput(taskMuonsHF, 0, mgr->GetCommonOutputContainer());

  char *fileName = "muonsHF.root";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosList",TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
  mgr->ConnectOutput(taskMuonsHF,1,coutput1);

  return taskMuonsHF;
}
