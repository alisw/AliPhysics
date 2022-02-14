AliAnalysisTask *AddTaskHFEemcQA(
                                 Bool_t UseTender=kTRUE,
                                 Bool_t FillElecSparse=kFALSE,
                                 Bool_t ClsTypeEMC=kTRUE, Bool_t ClsTypeDCAL=kTRUE,
                                 Bool_t hasTwoEMCTrigThres=kFALSE,
                                 Int_t  BitOption=0,
                                 Int_t MimCent = -1, Int_t MaxCent = -1,
                                 Int_t thEG1ADC=140, Int_t thEG2ADC=89,
                                 TString ContNameExt = "", TString centrality="V0M", Bool_t Wevt=kFALSE)
{
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
        if(MimCent==-1)
        {
            sprintf(calib,"wTender");
        }
        else
        {
            sprintf(calib,"wTender_%d_%d",MimCent,MaxCent);
        }
    }
    else
    {
        if(MimCent==-1)
        {
            sprintf(calib,"woTender");
        }
        else
        {
            sprintf(calib,"woTender_%d_%d",MimCent,MaxCent);
        }
    }
    
    if(ClsTypeEMC && !ClsTypeDCAL)ContNameExt+="_EMC";
    if(!ClsTypeEMC && ClsTypeDCAL)ContNameExt+="_DCAL";
    
    //Event plane task
    AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
    eventplaneTask->SelectCollisionCandidates(AliVEvent::kINT8 | AliVEvent::kINT7 | AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMCEGA | AliVEvent::kEMCEJE);
    if (type == "AOD") eventplaneTask->SetInput("AOD");
    eventplaneTask->SetTrackType("TPC");
    eventplaneTask->SetUsePtWeight();
    eventplaneTask->SetUsePhiWeight();
    eventplaneTask->SetSaveTrackContribution();
    mgr->AddTask(eventplaneTask);
    
    
    TString containerNameTPCEP = mgr->GetCommonFileName();
    containerNameTPCEP += ":PWGHF_hfeHFEemcQAEventPlane";
    containerNameTPCEP += ContNameExt;
    //TString SubcontainerNameTPCEP = Form("HFEemcQAINT8_%s",calib);
    TString SubcontainerNameTPCEP = Form("HFEemcQAINT8EventPlane_%s",calib);
    SubcontainerNameTPCEP += ContNameExt;
    
    AliAnalysisDataContainer *cinputTPCEP = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputTPCEP = mgr->CreateContainer(SubcontainerNameTPCEP,TList::Class(), AliAnalysisManager::kOutputContainer,
                                                                  containerNameTPCEP.Data());
    
    mgr->ConnectInput(eventplaneTask, 0,cinputTPCEP);
    mgr->ConnectOutput(eventplaneTask,1,coutputTPCEP);
    
    
    // +++ EMCal MB
    // INT7
    AliAnalysisTaskHFEemcQA *hfecalqa7 = new AliAnalysisTaskHFEemcQA("emcqa");
    mgr->AddTask(hfecalqa7);
    hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
    hfecalqa7->SetElecIDsparse(FillElecSparse);
    hfecalqa7->SetTenderSwitch(UseTender);
    hfecalqa7->SetThresholdEG1(thEG1ADC);
    hfecalqa7->SetThresholdEG2(thEG2ADC);
    hfecalqa7->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqa7->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqa7->SetCentralityMim(MimCent);
    hfecalqa7->SetCentralityMax(MaxCent);
    hfecalqa7->SetFilterBitOption(BitOption);
    hfecalqa7->SetCentralityEstimator(centrality.Data());
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PWGHF_hfeHFEemcQAINT7";
    containerName7 += ContNameExt;
    TString SubcontainerName7 = Form("HFEemcQAINT7_%s",calib);
    SubcontainerName7 += ContNameExt;
    AliAnalysisDataContainer *cinput7  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(SubcontainerName7, TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    mgr->ConnectInput(hfecalqa7, 0, cinput7);
    mgr->ConnectOutput(hfecalqa7, 1, coutput7);

    // EMCal EGA
    if(hasTwoEMCTrigThres)
    {
        // EMCal EGA EG1
        if(ClsTypeEMC){
        AliAnalysisTaskHFEemcQA *hfecalqaTrig01 = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfecalqaTrig01);
        hfecalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig01->SetEMCalTriggerEG1(kTRUE);
        hfecalqaTrig01->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig01->SetTenderSwitch(UseTender);
        hfecalqaTrig01->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig01->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig01->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaTrig01->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaTrig01->SetCentralityMim(MimCent);
        hfecalqaTrig01->SetCentralityMax(MaxCent);
        hfecalqaTrig01->SetFilterBitOption(BitOption);
        hfecalqaTrig01->SetCentralityEstimator(centrality.Data());
        
        TString containerName01 = mgr->GetCommonFileName();
        containerName01 += ":PWGHF_hfeHFEemcQATrigGAEG1";
        containerName01 += ContNameExt;
        TString SubcontainerName01 = Form("HFEemcQATrigGAEG1_%s",calib);
        SubcontainerName01 += ContNameExt;
        AliAnalysisDataContainer *cinput01  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput01 = mgr->CreateContainer(SubcontainerName01, TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
        mgr->ConnectInput(hfecalqaTrig01, 0, cinput01);
        mgr->ConnectOutput(hfecalqaTrig01, 1, coutput01);
        
        // EMCal EGA EG2
        AliAnalysisTaskHFEemcQA *hfecalqaTrig02 = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfecalqaTrig02);
        hfecalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig02->SetEMCalTriggerEG2(kTRUE);
        hfecalqaTrig02->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig02->SetTenderSwitch(UseTender);
        hfecalqaTrig02->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig02->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig02->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaTrig02->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaTrig02->SetCentralityMim(MimCent);
        hfecalqaTrig02->SetCentralityMax(MaxCent);
        hfecalqaTrig02->SetFilterBitOption(BitOption);
        hfecalqaTrig02->SetCentralityEstimator(centrality.Data());
        
        TString containerName02 = mgr->GetCommonFileName();
        containerName02 += ":PWGHF_hfeHFEemcQATrigGAEG2";
        containerName02 += ContNameExt;
        TString SubcontainerName02 = Form("HFEemcQATrigGAEG2_%s",calib);
        SubcontainerName02 += ContNameExt;
        AliAnalysisDataContainer *cinput02  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput02 = mgr->CreateContainer(SubcontainerName02, TList::Class(),AliAnalysisManager::kOutputContainer, containerName02.Data());
        mgr->ConnectInput(hfecalqaTrig02, 0, cinput02);
        mgr->ConnectOutput(hfecalqaTrig02, 1, coutput02);
        }
        if(ClsTypeDCAL){
        // DCal EGA DG1
        AliAnalysisTaskHFEemcQA *hfdcalqaTrig01 = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfdcalqaTrig01);
        hfdcalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfdcalqaTrig01->SetEMCalTriggerDG1(kTRUE);
        hfdcalqaTrig01->SetElecIDsparse(FillElecSparse);
        hfdcalqaTrig01->SetTenderSwitch(UseTender);
        hfdcalqaTrig01->SetThresholdEG1(thEG1ADC);
        hfdcalqaTrig01->SetThresholdEG2(thEG2ADC);
        hfdcalqaTrig01->SetClusterTypeEMC(ClsTypeEMC);
        hfdcalqaTrig01->SetClusterTypeDCAL(ClsTypeDCAL);
        hfdcalqaTrig01->SetCentralityMim(MimCent);
        hfdcalqaTrig01->SetCentralityMax(MaxCent);
        hfdcalqaTrig01->SetFilterBitOption(BitOption);
        hfdcalqaTrig01->SetCentralityEstimator(centrality.Data());
        
        TString containerName03 = mgr->GetCommonFileName();
        containerName03 += ":PWGHF_hfeHFEemcQATrigGADG1";
        containerName03 += ContNameExt;
        TString SubcontainerName03 = Form("HFEemcQATrigGADG1_%s",calib);
        SubcontainerName03 += ContNameExt;
        AliAnalysisDataContainer *cinput03  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput03 = mgr->CreateContainer(SubcontainerName03, TList::Class(),AliAnalysisManager::kOutputContainer, containerName03.Data());
        mgr->ConnectInput(hfdcalqaTrig01, 0, cinput03);
        mgr->ConnectOutput(hfdcalqaTrig01, 1, coutput03);
        
        // DCal EGA DG2
        AliAnalysisTaskHFEemcQA *hfdcalqaTrig02 = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfdcalqaTrig02);
        hfdcalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfdcalqaTrig02->SetEMCalTriggerDG2(kTRUE);
        hfdcalqaTrig02->SetElecIDsparse(FillElecSparse);
        hfdcalqaTrig02->SetTenderSwitch(UseTender);
        hfdcalqaTrig02->SetThresholdEG1(thEG1ADC);
        hfdcalqaTrig02->SetThresholdEG2(thEG2ADC);
        hfdcalqaTrig02->SetClusterTypeEMC(ClsTypeEMC);
        hfdcalqaTrig02->SetClusterTypeDCAL(ClsTypeDCAL);
        hfdcalqaTrig02->SetCentralityMim(MimCent);
        hfdcalqaTrig02->SetCentralityMax(MaxCent);
        hfdcalqaTrig02->SetFilterBitOption(BitOption);
        hfdcalqaTrig02->SetCentralityEstimator(centrality.Data());
        
        TString containerName04 = mgr->GetCommonFileName();
        containerName04 += ":PWGHF_hfeHFEemcQATrigGADG2";
        containerName04 += ContNameExt;
        TString SubcontainerName04 = Form("HFEemcQATrigGADG2_%s",calib);
        SubcontainerName04 += ContNameExt;
        AliAnalysisDataContainer *cinput04  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput04 = mgr->CreateContainer(SubcontainerName04, TList::Class(),AliAnalysisManager::kOutputContainer, containerName04.Data());
        mgr->ConnectInput(hfdcalqaTrig02, 0, cinput04);
        mgr->ConnectOutput(hfdcalqaTrig02, 1, coutput04);
        }
    }
    if(!hasTwoEMCTrigThres)
    {
        // EMCal EGA
        AliAnalysisTaskHFEemcQA *hfecalqaTrig0 = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfecalqaTrig0);
        hfecalqaTrig0->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig0->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig0->SetTenderSwitch(UseTender);
        hfecalqaTrig0->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig0->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig0->SetClusterTypeEMC(ClsTypeEMC);
        if(ClsTypeEMC)hfecalqaTrig0->SetEMCalTriggerEG1(kTRUE);
        hfecalqaTrig0->SetClusterTypeDCAL(ClsTypeDCAL);
        if(ClsTypeDCAL)hfecalqaTrig0->SetEMCalTriggerDG1(kTRUE);
        hfecalqaTrig0->SetCentralityMim(MimCent);
        hfecalqaTrig0->SetCentralityMax(MaxCent);
        hfecalqaTrig0->SetFilterBitOption(BitOption);
        hfecalqaTrig0->SetCentralityEstimator(centrality.Data());
        
        TString containerName1 = mgr->GetCommonFileName();
        containerName1 += ":PWGHF_hfeHFEemcQATrigGA";
        containerName1 += ContNameExt;
        TString SubcontainerName1 = Form("HFEemcQATrigGAEG_%s",calib);
        SubcontainerName1 += ContNameExt;
        AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName1, TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
        mgr->ConnectInput(hfecalqaTrig0, 0, cinput1);
        mgr->ConnectOutput(hfecalqaTrig0, 1, coutput1);
    }

        // central
        AliAnalysisTaskHFEemcQA *hfcalqaCentral = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfcalqaCentral);
        hfcalqaCentral->SelectCollisionCandidates(AliVEvent::kCentral);
        hfcalqaCentral->SetElecIDsparse(FillElecSparse);
        hfcalqaCentral->SetTenderSwitch(UseTender);
        hfcalqaCentral->SetThresholdEG1(thEG1ADC);
        hfcalqaCentral->SetThresholdEG2(thEG2ADC);
        hfcalqaCentral->SetClusterTypeEMC(ClsTypeEMC);
        hfcalqaCentral->SetClusterTypeDCAL(ClsTypeDCAL);
        hfcalqaCentral->SetCentralityMim(MimCent);
        hfcalqaCentral->SetCentralityMax(MaxCent);
        hfcalqaCentral->SetFilterBitOption(BitOption);
        hfcalqaCentral->SetCentralityEstimator(centrality.Data());
        
        TString containerName05 = mgr->GetCommonFileName();
        containerName05 += ":PWGHF_hfeHFEemcQACentral";
        containerName05 += ContNameExt;
        TString SubcontainerName05 = Form("HFEemcQAcentral_%s",calib);
        SubcontainerName05 += ContNameExt;
        AliAnalysisDataContainer *cinput05  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput05 = mgr->CreateContainer(SubcontainerName05, TList::Class(),AliAnalysisManager::kOutputContainer, containerName05.Data());
        mgr->ConnectInput(hfcalqaCentral, 0, cinput05);
        mgr->ConnectOutput(hfcalqaCentral, 1, coutput05);


        // semi-central
        AliAnalysisTaskHFEemcQA *hfcalqaSemiCentral = new AliAnalysisTaskHFEemcQA("emcqa");
        mgr->AddTask(hfcalqaSemiCentral);
        hfcalqaSemiCentral->SelectCollisionCandidates(AliVEvent::kSemiCentral);
        hfcalqaSemiCentral->SetElecIDsparse(FillElecSparse);
        hfcalqaSemiCentral->SetTenderSwitch(UseTender);
        hfcalqaSemiCentral->SetThresholdEG1(thEG1ADC);
        hfcalqaSemiCentral->SetThresholdEG2(thEG2ADC);
        hfcalqaSemiCentral->SetClusterTypeEMC(ClsTypeEMC);
        hfcalqaSemiCentral->SetClusterTypeDCAL(ClsTypeDCAL);
        hfcalqaSemiCentral->SetCentralityMim(MimCent);
        hfcalqaSemiCentral->SetCentralityMax(MaxCent);
        hfcalqaSemiCentral->SetFilterBitOption(BitOption);
        hfcalqaSemiCentral->SetCentralityEstimator(centrality.Data());
        
        TString containerName06 = mgr->GetCommonFileName();
        containerName06 += ":PWGHF_hfeHFEemcQASemiCentral";
        containerName06 += ContNameExt;
        TString SubcontainerName06 = Form("HFEemcQASemicentral_%s",calib);
        SubcontainerName06 += ContNameExt;
        AliAnalysisDataContainer *cinput06  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput06 = mgr->CreateContainer(SubcontainerName06, TList::Class(),AliAnalysisManager::kOutputContainer, containerName06.Data());
        mgr->ConnectInput(hfcalqaSemiCentral, 0, cinput06);
        mgr->ConnectOutput(hfcalqaSemiCentral, 1, coutput06);

    if(Wevt)
      {
       AliAnalysisTaskHFEemcQA *hfecalqaWmc = new AliAnalysisTaskHFEemcQA("emcqa");
       mgr->AddTask(hfecalqaWmc);
       hfecalqaWmc->SelectCollisionCandidates(AliVEvent::kINT7);
       hfecalqaWmc->SetElecIDsparse(FillElecSparse);
       hfecalqaWmc->SetTenderSwitch(UseTender);
       hfecalqaWmc->SetThresholdEG1(thEG1ADC);
       hfecalqaWmc->SetThresholdEG2(thEG2ADC);
       hfecalqaWmc->SetClusterTypeEMC(ClsTypeEMC);
       hfecalqaWmc->SetClusterTypeDCAL(ClsTypeDCAL);
       hfecalqaWmc->SetCentralityMim(MimCent);
       hfecalqaWmc->SetCentralityMax(MaxCent);
       hfecalqaWmc->SetFilterBitOption(BitOption);
       hfecalqaWmc->SetCentralityEstimator(centrality.Data());
       hfecalqaWmc->SetMC_W(Wevt);
    
       TString containerNameW = mgr->GetCommonFileName();
       containerNameW += ":PWGHF_hfeHFEemcQAWmc";
       containerNameW += ContNameExt;
       TString SubcontainerNameW = Form("HFEemcQAWmc_%s",calib);
       SubcontainerNameW += ContNameExt;
       AliAnalysisDataContainer *cinputW  = mgr->GetCommonInputContainer();
       AliAnalysisDataContainer *coutputW = mgr->CreateContainer(SubcontainerNameW, TList::Class(),AliAnalysisManager::kOutputContainer, containerNameW.Data());
       mgr->ConnectInput(hfecalqaWmc, 0, cinputW);
       mgr->ConnectOutput(hfecalqaWmc, 1, coutputW);
     }

  /*
    // EMCal EJE
    AliAnalysisTaskHFEemcQA *hfecalqaTrig1 = new AliAnalysisTaskHFEemcQA("emcqa");
    mgr->AddTask(hfecalqaTrig1);
    hfecalqaTrig1->SelectCollisionCandidates(AliVEvent::kEMCEJE);
    hfecalqaTrig1->SetElecIDsparse(FillElecSparse);
    hfecalqaTrig1->SetTenderSwitch(UseTender);
    hfecalqaTrig1->SetThresholdEG1(thEG1ADC);
    hfecalqaTrig1->SetThresholdEG2(thEG2ADC);
    hfecalqaTrig1->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqaTrig1->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqaTrig1->SetCentralityMim(MimCent);
    hfecalqaTrig1->SetCentralityMax(MaxCent);
    hfecalqaTrig1->SetCentralityEstimator(centrality.Data());
    
    TString containerName2 = mgr->GetCommonFileName();
    containerName2 += ":PWGHF_hfeHFEemcQATrigJE";
    containerName2 += ContNameExt;
    TString SubcontainerName2 = Form("HFEemcQATrigJE_%s",calib);
    SubcontainerName2 += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName2, TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
    mgr->ConnectInput(hfecalqaTrig1, 0, cinput);
    mgr->ConnectOutput(hfecalqaTrig1, 1, coutput1);
 */   
    //return hfecalqa;

    return NULL;
}
