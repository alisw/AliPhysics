AliAnalysisTask *AddTaskFlowTPCEMCalEP( 
 Double_t fP2_lowPtEta0010 = 1., Double_t fP3_lowPtEta0010 = 1., Double_t fP4_lowPtEta0010 = 1., 
 Double_t fP2_highPtEta0010 = 1., Double_t fP3_highPtEta0010 = 1., Double_t fP4_highPtEta0010 = 1., 
 Double_t fP2_lowPtPi00010 = 1., Double_t fP3_lowPtPi00010 = 1., Double_t fP4_lowPtPi00010 = 1., 
 Double_t fP2_highPtPi00010 = 1., Double_t fP3_highPtPi00010 = 1., Double_t fP4_highPtPi00010 = 1., 
 Double_t fP2_lowPtEta1020 = 1., Double_t fP3_lowPtEta1020 = 1., Double_t fP4_lowPtEta1020 = 1., 
 Double_t fP2_highPtEta1020 = 1., Double_t fP3_highPtEta1020 = 1., Double_t fP4_highPtEta1020 = 1., 
 Double_t fP2_lowPtPi01020 = 1., Double_t fP3_lowPtPi01020 = 1., Double_t fP4_lowPtPi01020 = 1., 
 Double_t fP2_highPtPi01020 = 1., Double_t fP3_highPtPi01020 = 1., Double_t fP4_highPtPi01020 = 1., 
 Double_t fP2_lowPtEta2040 = 1., Double_t fP3_lowPtEta2040 = 1., Double_t fP4_lowPtEta2040 = 1., 
 Double_t fP2_highPtEta2040 = 1., Double_t fP3_highPtEta2040 = 1., Double_t fP4_highPtEta2040 = 1., 
 Double_t fP2_lowPtPi02040 = 1., Double_t fP3_lowPtPi02040 = 1., Double_t fP4_lowPtPi02040 = 1., 
 Double_t fP2_highPtPi02040 = 1., Double_t fP3_highPtPi02040 = 1., Double_t fP4_highPtPi02040 = 1.,
 TString ID="ContName")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowTPCEMCalEP", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowTPCEMCalEP", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskFlowTPCEMCalEP", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  
  //Event plane task
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kEMCEGA | AliVEvent::kEMCEJE);

  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();

  mgr->AddTask(eventplaneTask);

  TString containerName0 = mgr->GetCommonFileName();
  containerName0 += ":PWGHF_hfeCalEventPlane";
  containerName0 += ID;
  
  TString name0 = "EPStat";
  name0 += ID;
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name0.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,containerName0.Data());
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);

  //analysis task 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TPCEMCal_EP.C");

  AliAnalysisTaskFlowTPCEMCalEP *taskMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,fP2_lowPtEta0010,fP3_lowPtEta0010,fP4_lowPtEta0010,fP2_highPtEta0010,fP3_highPtEta0010,fP4_highPtEta0010,fP2_lowPtPi00010,fP3_lowPtPi00010,fP4_lowPtPi00010,fP2_highPtPi00010,fP3_highPtPi00010,fP4_highPtPi00010,fP2_lowPtEta1020,fP3_lowPtEta1020,fP4_lowPtEta1020,fP2_highPtEta1020,fP3_highPtEta1020,fP4_highPtEta1020,fP2_lowPtPi01020,fP3_lowPtPi01020,fP4_lowPtPi01020,fP2_highPtPi01020,fP3_highPtPi01020,fP4_highPtPi01020,fP2_lowPtEta2040,fP3_lowPtEta2040,fP4_lowPtEta2040,fP2_highPtEta2040,fP3_highPtEta2040,fP4_highPtEta2040,fP2_lowPtPi02040,fP3_lowPtPi02040,fP4_lowPtPi02040,fP2_highPtPi02040,fP3_highPtPi02040,fP4_highPtPi02040);
  AliAnalysisTaskFlowTPCEMCalEP *taskcorrMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,fP2_lowPtEta0010,fP3_lowPtEta0010,fP4_lowPtEta0010,fP2_highPtEta0010,fP3_highPtEta0010,fP4_highPtEta0010,fP2_lowPtPi00010,fP3_lowPtPi00010,fP4_lowPtPi00010,fP2_highPtPi00010,fP3_highPtPi00010,fP4_highPtPi00010,fP2_lowPtEta1020,fP3_lowPtEta1020,fP4_lowPtEta1020,fP2_highPtEta1020,fP3_highPtEta1020,fP4_highPtEta1020,fP2_lowPtPi01020,fP3_lowPtPi01020,fP4_lowPtPi01020,fP2_highPtPi01020,fP3_highPtPi01020,fP4_highPtPi01020,fP2_lowPtEta2040,fP3_lowPtEta2040,fP4_lowPtEta2040,fP2_highPtEta2040,fP3_highPtEta2040,fP4_highPtEta2040,fP2_lowPtPi02040,fP3_lowPtPi02040,fP4_lowPtPi02040,fP2_highPtPi02040,fP3_highPtPi02040,fP4_highPtPi02040);
  AliAnalysisTaskFlowTPCEMCalEP *taskTR = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,fP2_lowPtEta0010,fP3_lowPtEta0010,fP4_lowPtEta0010,fP2_highPtEta0010,fP3_highPtEta0010,fP4_highPtEta0010,fP2_lowPtPi00010,fP3_lowPtPi00010,fP4_lowPtPi00010,fP2_highPtPi00010,fP3_highPtPi00010,fP4_highPtPi00010,fP2_lowPtEta1020,fP3_lowPtEta1020,fP4_lowPtEta1020,fP2_highPtEta1020,fP3_highPtEta1020,fP4_highPtEta1020,fP2_lowPtPi01020,fP3_lowPtPi01020,fP4_lowPtPi01020,fP2_highPtPi01020,fP3_highPtPi01020,fP4_highPtPi01020,fP2_lowPtEta2040,fP3_lowPtEta2040,fP4_lowPtEta2040,fP2_highPtEta2040,fP3_highPtEta2040,fP4_highPtEta2040,fP2_lowPtPi02040,fP3_lowPtPi02040,fP4_lowPtPi02040,fP2_highPtPi02040,fP3_highPtPi02040,fP4_highPtPi02040);
 
  mgr->AddTask(taskcorrMB);
  mgr->AddTask(taskMB);
  mgr->AddTask(taskTR);
  
  // Flattened semi central trigger

  taskcorrMB->SelectCollisionCandidates(AliVEvent::kAny);

  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_hfeCalcorrSemiCentralV2";
  containerName1 += ID;
  
  TString name1 = "histcorrMB";
  name1 += ID;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  mgr->ConnectInput(taskcorrMB, 0, cinput);
  mgr->ConnectOutput(taskcorrMB, 1, coutput1);

  // Central trigger
  taskMB->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral);

  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalCentralV2";
  containerName2 += ID;
  
  TString name2 = "histMB";
  name2 += ID;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);
  
  //L1 gamma trigger
  taskTR->SelectCollisionCandidates(AliVEvent::kEMCEGA);

  TString containerName3 = mgr->GetCommonFileName();
  containerName3 += ":PWGHF_hfeCalL1GammaV2";
  containerName3 += ID;
  
  TString name3 = "histTR";
  name3 += ID;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name3.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName3.Data());
  mgr->ConnectInput(taskTR, 0, cinput);
  mgr->ConnectOutput(taskTR, 1, coutput1);
  
  if(MCthere){
    
    AliAnalysisTaskFlowTPCEMCalEP *taskMC = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,fP2_lowPtEta0010,fP3_lowPtEta0010,fP4_lowPtEta0010,fP2_highPtEta0010,fP3_highPtEta0010,fP4_highPtEta0010,fP2_lowPtPi00010,fP3_lowPtPi00010,fP4_lowPtPi00010,fP2_highPtPi00010,fP3_highPtPi00010,fP4_highPtPi00010,fP2_lowPtEta1020,fP3_lowPtEta1020,fP4_lowPtEta1020,fP2_highPtEta1020,fP3_highPtEta1020,fP4_highPtEta1020,fP2_lowPtPi01020,fP3_lowPtPi01020,fP4_lowPtPi01020,fP2_highPtPi01020,fP3_highPtPi01020,fP4_highPtPi01020,fP2_lowPtEta2040,fP3_lowPtEta2040,fP4_lowPtEta2040,fP2_highPtEta2040,fP3_highPtEta2040,fP4_highPtEta2040,fP2_lowPtPi02040,fP3_lowPtPi02040,fP4_lowPtPi02040,fP2_highPtPi02040,fP3_highPtPi02040,fP4_highPtPi02040);
    mgr->AddTask(taskMC);
    
    taskMC->SelectCollisionCandidates(AliVEvent::kMB);
    
    TString containerName4 = mgr->GetCommonFileName();
    containerName4 += ":PWGHF_hfeCalMCV2";
    containerName4 += ID;
    
    TString name4 = "histMC";
    name4 += ID;
  
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name4.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName4.Data());
    mgr->ConnectInput(taskMC, 0, cinput);
    mgr->ConnectOutput(taskMC, 1, coutput1);
  }
  
  
  return NULL;
}

