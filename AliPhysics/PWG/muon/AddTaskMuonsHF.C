AliAnalysisTaskSEMuonsHF* AddTaskMuonsHF(Int_t mode=0, Int_t passN=2, Bool_t isMC=kFALSE, Bool_t isTree=kFALSE, Bool_t isSel=kTRUE)
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
  Double_t cutsEvsH[5] ={-999999.0,   // low limit of Ncontrs
                          999999.0,   // up limit of |vz|
                          999999.0,   // up limit of vt
                         -999999.0,   // centrality min
                          999999.0};  // centrality max
  AliMuonTrackCuts *cutsMuon = new AliMuonTrackCuts("cutsMuon", "cutsMuon"); cutsMuon->SetIsMC(isMC); cutsMuon->SetPassNumber(passN);
  AliMuonPairCuts  *cutsDimu = new  AliMuonPairCuts("cutsDimu", "cutsDimu"); cutsDimu->SetIsMC(isMC); cutsDimu->GetMuonTrackCuts()->SetPassNumber(passN);
  AliAnalysisTaskSEMuonsHF *taskMuonsHF = new AliAnalysisTaskSEMuonsHF("MuonsHF Analysis Task",*cutsMuon,*cutsDimu);
  taskMuonsHF->SetAnaMode(mode);
  taskMuonsHF->SetUseMC(isMC);
  taskMuonsHF->SetIsOutputTree(isTree);
  taskMuonsHF->SetEvsHCuts(cutsEvsH);
  if (isSel) taskMuonsHF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);
  mgr->AddTask(taskMuonsHF);

  mgr->ConnectInput(taskMuonsHF, 0, mgr->GetCommonInputContainer());
  if (isTree) mgr->ConnectOutput(taskMuonsHF, 0, mgr->GetCommonOutputContainer());

  char *fileName = "muonsHF.root";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosList",TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
  mgr->ConnectOutput(taskMuonsHF,1,coutput1);

  return taskMuonsHF;
}
