AliAnalysisTask *AddTaskHFECalSys(int sysID, int TPCclust, int TPCclustPID, int Nits, int ITSstat, int QAhist)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEECal", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFECal", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFECal", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  
  //analysis task 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFECalSys.C");

  AliAnalysisTaskHFECal *hfetaskCent = ConfigHFECalSys(MCthere,TPCclust,TPCclustPID,Nits,ITSstat,QAhist);
  AliAnalysisTaskHFECal *hfetaskTrig= ConfigHFECalSys(MCthere,TPCclust,TPCclustPID,Nits,ITSstat,QAhist);
 
  mgr->AddTask(hfetaskCent);
  mgr->AddTask(hfetaskTrig);
  
  // Semi-central trigger
  hfetaskCent->SelectCollisionCandidates(AliVEvent::kCentral);
  
  TString Output1(TString::Format("HFE_Results_EMCALCentral%d",sysID));
  TString containerName = mgr->GetCommonFileName();
  TString appendix(TString::Format(":PWGHF_hfeCalCentral%d",sysID));
  containerName += appendix;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Output1, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(hfetaskCent, 0, cinput);
  mgr->ConnectOutput(hfetaskCent, 1, coutput1);
  
  //L1 gamma trigger
  hfetaskTrig->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  
  TString Output2(TString::Format("HFE_Results_EMCalTrigEGA%d",sysID));
  TString containerName2 = mgr->GetCommonFileName();
  TString appendix(TString::Format(":PWGHF_hfeCalTrigEGA%d",sysID));
  containerName2 += appendix;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Output2, TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(hfetaskTrig, 0, cinput);
  mgr->ConnectOutput(hfetaskTrig, 1, coutput1);
  
    //MB trigger
    AliAnalysisTaskHFECal *hfetaskMB = ConfigHFECalSys(MCthere,TPCclust,TPCclustPID,Nits,ITSstat,QAhist);
    mgr->AddTask(hfetaskMB);
    hfetaskMB->SelectCollisionCandidates(AliVEvent::kMB);

    TString Output3(TString::Format("HFE_Results_EMCalMB%d",sysID));
    TString containerName3 = mgr->GetCommonFileName();
    TString appendix(TString::Format(":PWGHF_hfeCalkMB%d",sysID));
    containerName3 += appendix;

     AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
     AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Output3, TList::Class(),AliAnalysisManager::kOutputContainer, containerName3.Data());
     mgr->ConnectInput(hfetaskMB, 0, cinput);
     mgr->ConnectOutput(hfetaskMB, 1, coutput1);
  

  return NULL;
}
