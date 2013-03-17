// $Id: AddTaskEMCALPi0V2.C 56081 2012-05-01 08:57:08Z loizides $

AliAnalysisTask *AddTaskEMCALPi0V2 (
  TString trackName  = "PicoTracks",
  Double_t Ecut      = 1,   
  Double_t M02cut    = 0.5, 
  Double_t fDrCut    = 0.025, 
  Bool_t IsV1cus     = 0,
  TString V1ClusName = "CaloClusters", 
  TString V2ClusName = "caloClusters", 
  TString trigClass  = "",
  Bool_t IsPhosCali  = kFALSE,
  Int_t EvtType      = 5 
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEMCALPi0V2", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALPi0V2", "This task requires an input event handler");
    return NULL;
  }

  TString Input;
  AliAnalysisTaskPi0V2* taskMB = new  AliAnalysisTaskPi0V2("Pi0v2Task");
  if(EvtType == 1){ //central
    taskMB->SelectCollisionCandidates(AliVEvent::kCentral);
    Input = "kCentral";
  } else if (EvtType == 2){ //SemiCentral
    taskMB->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    Input = "kSemiCentral";
  } else if (EvtType == 3){ //kMB 
    taskMB->SelectCollisionCandidates(AliVEvent::kMB);
    Input = "kMB";
  } else if (EvtType == 4){ //Central + SemiCentral 
    taskMB->SelectCollisionCandidates(AliVEvent::kCentral | AliVEvent::kSemiCentral);
    Input = "Central_SemiCentral";
  } else if (EvtType == 5){ //Central + SemiCentral + kMB
    taskMB->SelectCollisionCandidates(AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB);
    Input = "ALLMB";
  }
  taskMB->SetTracksName(trackName.Data());
  taskMB->SetClusE(Ecut);
  taskMB->SetClusM02(M02cut);
  taskMB->SetDrCut(fDrCut);
  taskMB->SetIsV1Clus(IsV1cus);
  taskMB->SetV1ClusName(V1ClusName);
  taskMB->SetV2ClusName(V2ClusName);
  taskMB->SetTrigClass(trigClass);
  taskMB->SetIsPHOSCali(IsPhosCali);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGGA_EMCalpi0v2";

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(
    Form("%s_E%1.2f_M02%1.2f", Input.Data(), Ecut, M02cut), 
    TList::Class(),
    AliAnalysisManager::kOutputContainer, 
    containerName.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput2);

  return taskMB;
}
