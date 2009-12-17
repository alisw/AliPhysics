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

  /*if (isMC) {
    AliMCEventHandler *mcH = mgr->GetMCtruthEventHandler();
    if (!mcH) {
      ::Error("AddTaskMuonsHF", "MuonsHF task needs the manager to have an MC evnet handler.");
      return NULL;
    }
  }*/

  if (isTree) {
    AliAODHandler *aodH = (AliAODHandler*)mgr->GetOutputEventHandler();
    if (!aodH) {
      ::Error("AddTaskMuonsHF", "MuonsHF task needs the manager to have an AOD output handler.");
      return NULL;
    }
  }

  // set cuts for single muon track selection
  Double_t cuts[10]={-1.,      // 0, min of 3-momentum
                     999999.,  // 1, max of 3-momnentum
                     1.,       // 2, PtMin
                     999999.,  // 3, PtMax
                     -4.,      // 4, EtaMin
                     -2.5,     // 5, EtaMax
                     -1.,      // 6, DCAmin
                     10.,      // 7, DCAmax
                     0.5,      // 8, for trigger
                     3.5.      // 9, for trigger
                    };

  AliAnalysisTaskSEMuonsHF *taskMuonsHF = new AliAnalysisTaskSEMuonsHF("MuonsHF Analysis Task");
  taskMuonsHF->SetAnaMode(mode);
  taskMuonsHF->SetIsUseMC(isMC);
  taskMuonsHF->SetIsOutputTree(isTree);
  taskMuonsHF->SetSingleMuonCuts(cuts);
  mgr->AddTask(taskMuonsHF);
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3Muon_MuonHF";

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("listHisEventH",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("listHisSingleMuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("listHisDimuon",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  mgr->ConnectInput(taskMuonsHF,0,mgr->GetCommonInputContainer());
  if (isTree) mgr->ConnectOutput(taskMuonsHF,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(taskMuonsHF,1,coutput2);
  mgr->ConnectOutput(taskMuonsHF,2,coutput3);
  mgr->ConnectOutput(taskMuonsHF,3,coutput4);

  return taskMuonsHF;
}
