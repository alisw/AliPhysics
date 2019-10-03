AliAnalysisTask *AddTaskEMCALPi0V2hardCodeEP(Double_t EvtMthod=2, TString trackName="PicoTrack", 
					     Double_t Ecut = 1,   Double_t M02cut = 0.5, Double_t fDrCut=0.025, Bool_t IsV1cus = 0,
					     TString V1ClusName="CaloCluster", TString V2ClusName="CaloCluster"
				            )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEMCALPi0V2hardCodeEP", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALPi0V2hardCodeEP", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskEMCALPi0V2hardCodeEP", "The tasks exits because AODs are in input");
    return NULL;
  }
  
  //Event plane task
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();
  
  AliESDtrackCuts* epTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
  epTrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  epTrackCuts->SetMinNClustersTPC(50);
  epTrackCuts->SetMaxChi2PerClusterTPC(4);
  epTrackCuts->SetAcceptKinkDaughters(kFALSE);
  epTrackCuts->SetRequireTPCRefit(kTRUE);
  epTrackCuts->SetMaxDCAToVertexZ(3.2);
  epTrackCuts->SetMaxDCAToVertexXY(2.4);
  epTrackCuts->SetPtRange(0.15, 20);
  eventplaneTask->SetPersonalESDtrackCuts(epTrackCuts);

  mgr->AddTask(eventplaneTask);

  TString containerName3 = mgr->GetCommonFileName();
  containerName3 += ":PWGGA_pi0v2CalEventPlane";
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("EPStatTPC_E%1.2f_M02%1.2f", Ecut, M02cut),TList::Class(), AliAnalysisManager::kOutputContainer,containerName3.Data());
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);

  //analysis task 

  AliAnalysisTaskPi0V2* taskMB = new  AliAnalysisTaskPi0V2("Pi0v2Task"); 
  taskMB->SetEventMethod(EvtMthod);
  taskMB->SetTracksName(trackName.Data());
  taskMB->SetClusE(Ecut);
  taskMB->SetClusM02(M02cut);
  taskMB->SetDrCut(fDrCut);
  taskMB->SetIsV1Clus(IsV1cus);
  taskMB->SetV1ClusName(V1ClusName);
  taskMB->SetV2ClusName(V2ClusName);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGGA_pi0v2CalSemiCentral";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("histv2task_E%1.2f_M02%1.2f", Ecut, M02cut), TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput2);
  
  return NULL;
}
