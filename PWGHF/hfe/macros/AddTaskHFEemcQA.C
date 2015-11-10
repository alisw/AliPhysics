AliAnalysisTask *AddTaskHFEemcQA(Bool_t UseTender=kTRUE, Bool_t FillElecSparse=kFALSE, Bool_t ispPb=kFALSE,Int_t thEG1ADC=140,Int_t thEG2ADC=89){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFE", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFE", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  char calib[100];
  if(UseTender)
    {
      sprintf(calib,"wTender");
    }
  else
    {
      sprintf(calib,"woTender");
    }
  
  // +++ EMCal MB
  AliAnalysisTaskHFEemcQA *hfecalqa = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqa);
  hfecalqa->SelectCollisionCandidates(AliVEvent::kINT8);
  hfecalqa->SetElecIDsparse(FillElecSparse);
  hfecalqa->SetTenderSwitch(UseTender);
  hfecalqa->SetEMCalTriggerEG1(thEG1ADC);
  hfecalqa->SetEMCalTriggerEG2(thEG2ADC);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeHFEemcQAINT8";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQAINT8_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(hfecalqa, 0, cinput);
  mgr->ConnectOutput(hfecalqa, 1, coutput1); 
  
  AliAnalysisTaskHFEemcQA *hfecalqa7 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqa7);
  hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
  hfecalqa7->SetElecIDsparse(FillElecSparse);
  hfecalqa7->SetTenderSwitch(UseTender);
  hfecalqa7->SetEMCalTriggerEG1(thEG1ADC);
  hfecalqa7->SetEMCalTriggerEG2(thEG2ADC);
  
  TString containerName7 = mgr->GetCommonFileName();
  containerName7 += ":PWGHF_hfeHFEemcQAINT7";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQAINT7_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
  mgr->ConnectInput(hfecalqa7, 0, cinput);
  mgr->ConnectOutput(hfecalqa7, 1, coutput1); 
 
  // EMCal L0
  // + kEMC7
  AliAnalysisTaskHFEemcQA *hfecalqaL07 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqaL07);
  hfecalqaL07->SelectCollisionCandidates(AliVEvent::kEMC7);
  hfecalqaL07->SetElecIDsparse(FillElecSparse);
  hfecalqaL07->SetTenderSwitch(UseTender);
  hfecalqaL07->SetEMCalTriggerEG1(thEG1ADC);
  hfecalqaL07->SetEMCalTriggerEG2(thEG2ADC);

  TString containerNameL07 = mgr->GetCommonFileName();
  containerNameL07 += ":PWGHF_hfeHFEemcQAEMC7";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQAEMC7_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerNameL07.Data());
  mgr->ConnectInput(hfecalqaL07, 0, cinput);
  mgr->ConnectOutput(hfecalqaL07, 1, coutput1); 

  // + kEMC8
  AliAnalysisTaskHFEemcQA *hfecalqaL08 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqaL08);
  hfecalqaL08->SelectCollisionCandidates(AliVEvent::kEMC8);
  hfecalqaL08->SetElecIDsparse(FillElecSparse);
  hfecalqaL08->SetTenderSwitch(UseTender);
  hfecalqaL08->SetEMCalTriggerEG1(thEG1ADC);
  hfecalqaL08->SetEMCalTriggerEG2(thEG2ADC);
  
  TString containerNameL08 = mgr->GetCommonFileName();
  containerNameL08 += ":PWGHF_hfeHFEemcQAEMC8";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQAEMC8_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerNameL08.Data());
  mgr->ConnectInput(hfecalqaL08, 0, cinput);
  mgr->ConnectOutput(hfecalqaL08, 1, coutput1); 
 

  // EMCal EGA
  if(ispPb)
    {
      // EMCal EGA EG1 
      AliAnalysisTaskHFEemcQA *hfecalqaTrig01 = new AliAnalysisTaskHFEemcQA("emcqa");
      mgr->AddTask(hfecalqaTrig01);
      hfecalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      hfecalqaTrig01->SetEMCalTriggerEG1(kTRUE);
      hfecalqaTrig01->SetElecIDsparse(FillElecSparse);
      hfecalqaTrig01->SetTenderSwitch(UseTender);
      hfecalqaTrig01->SetEMCalTriggerEG1(thEG1ADC);
      hfecalqaTrig01->SetEMCalTriggerEG2(thEG2ADC);

      TString containerName01 = mgr->GetCommonFileName();
      containerName01 += ":PWGHF_hfeHFEemcQATrigGAEG1";
      AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQATrigGAEG1_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
      mgr->ConnectInput(hfecalqaTrig01, 0, cinput);
      mgr->ConnectOutput(hfecalqaTrig01, 1, coutput1); 
      
      // EMCal EGA EG2 
      AliAnalysisTaskHFEemcQA *hfecalqaTrig02 = new AliAnalysisTaskHFEemcQA("emcqa");
      mgr->AddTask(hfecalqaTrig02);
      hfecalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      hfecalqaTrig02->SetEMCalTriggerEG2(kTRUE);
      hfecalqaTrig02->SetElecIDsparse(FillElecSparse);
      hfecalqaTrig02->SetTenderSwitch(UseTender);
      hfecalqaTrig02->SetEMCalTriggerEG1(thEG1ADC);
	  hfecalqaTrig02->SetEMCalTriggerEG2(thEG2ADC);
  
      TString containerName02 = mgr->GetCommonFileName();
      containerName02 += ":PWGHF_hfeHFEemcQATrigGAEG2";
      AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQATrigGAEG2_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName02.Data());
      mgr->ConnectInput(hfecalqaTrig02, 0, cinput);
      mgr->ConnectOutput(hfecalqaTrig02, 1, coutput1); 
    }
  if(!ispPb)
    {
      // EMCal EGA
      AliAnalysisTaskHFEemcQA *hfecalqaTrig0 = new AliAnalysisTaskHFEemcQA("emcqa");
      mgr->AddTask(hfecalqaTrig0);
      hfecalqaTrig0->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      hfecalqaTrig0->SetElecIDsparse(FillElecSparse);
      hfecalqaTrig0->SetTenderSwitch(UseTender);
      hfecalqaTrig0->SetEMCalTriggerEG1(thEG1ADC);
      hfecalqaTrig0->SetEMCalTriggerEG2(thEG2ADC);

      TString containerName1 = mgr->GetCommonFileName();
      containerName1 += ":PWGHF_hfeHFEemcQATrigGA";
      AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQATrigGA_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
      mgr->ConnectInput(hfecalqaTrig0, 0, cinput);
      mgr->ConnectOutput(hfecalqaTrig0, 1, coutput1); 
    }
  
  // EMCal EJE
  AliAnalysisTaskHFEemcQA *hfecalqaTrig1 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqaTrig1);
  hfecalqaTrig1->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  hfecalqaTrig1->SetElecIDsparse(FillElecSparse);
  hfecalqaTrig1->SetTenderSwitch(UseTender);
  hfecalqaTrig1->SetEMCalTriggerEG1(thEG1ADC);
  hfecalqaTrig1->SetEMCalTriggerEG2(thEG2ADC);
  
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeHFEemcQATrigJE";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFEemcQATrigJE_%s",calib), TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(hfecalqaTrig1, 0, cinput);
  mgr->ConnectOutput(hfecalqaTrig1, 1, coutput1); 
  
  //return hfecalqa;
  return NULL;
}
