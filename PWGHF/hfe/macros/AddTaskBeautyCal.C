AliAnalysisTask *AddTaskBeautyCal(
                                 Bool_t UseTender=kTRUE,
                                 Bool_t FillElecSparse=kFALSE,
                                 Bool_t ClsTypeEMC=kTRUE, Bool_t ClsTypeDCAL=kTRUE,
                                 Bool_t FillElec=kFALSE,
                                 Bool_t ispPb=kFALSE,
                                 Int_t MimCent = -1, Int_t MaxCent = -1, 
                                 Double_t nSig = -1.0,
                                 Int_t EtaRange = 0,
                                 Int_t ITSchi2 = 26,
                                 Bool_t pileup = kFALSE,
                                 Double_t mass = 0.1, Double_t ptAss = 0.2, 
                                 Int_t thEG1ADC=140, Int_t thEG2ADC=89,
                                 TString ContNameExt = "")
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

    // +++ EMCal MB
    // INT8
    /*
    AliAnalysisTaskBeautyCal *hfecalqa = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqa);
    hfecalqa->SelectCollisionCandidates(AliVEvent::kINT8);
    hfecalqa->SetElecIDsparse(FillElecSparse);
    hfecalqa->SetElecSys(FillElec);
    hfecalqa->SetTenderSwitch(UseTender);
    hfecalqa->SetThresholdEG1(thEG1ADC);
    hfecalqa->SetThresholdEG2(thEG2ADC);
    hfecalqa->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqa->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqa->SetCentralityMim(MimCent);
    hfecalqa->SetCentralityMax(MaxCent);
    hfecalqa->SetInvMassCut0(mass);
    hfecalqa->SetInvMassCut1(ptAss);
    hfecalqa->SetEtaRange(EtaRange);
    
    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGHF_hfeBeautyCalINT8";
    containerName += ContNameExt;
    TString SubcontainerName = Form("BeautyCalINT8_%s",calib);
    SubcontainerName += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(hfecalqa, 0, cinput);
    mgr->ConnectOutput(hfecalqa, 1, coutput1);
    */
    // INT7
    AliAnalysisTaskBeautyCal *hfecalqa7 = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqa7);
    hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
    hfecalqa7->SetElecIDsparse(FillElecSparse);
    hfecalqa7->SetElecSys(FillElec);
    hfecalqa7->SetTenderSwitch(UseTender);
    hfecalqa7->SetThresholdEG1(thEG1ADC);
    hfecalqa7->SetThresholdEG2(thEG2ADC);
    hfecalqa7->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqa7->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqa7->SetCentralityMim(MimCent);
    hfecalqa7->SetCentralityMax(MaxCent);
    hfecalqa7->SetInvMassCut0(mass);
    hfecalqa7->SetInvMassCut1(ptAss);
    hfecalqa7->SetEtaRange(EtaRange);
    hfecalqa7->SetITSchi2(ITSchi2);
    hfecalqa7->SetMinSig(nSig);
    hfecalqa7->SetPileUpCut(pileup);
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PWGHF_hfeBeautyCalINT7";
    containerName7 += ContNameExt;
    TString SubcontainerName7 = Form("BeautyCalINT7_%s",calib);
    SubcontainerName7 += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName7, TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    mgr->ConnectInput(hfecalqa7, 0, cinput);
    mgr->ConnectOutput(hfecalqa7, 1, coutput1);
    
    // EMCal L0
    // + kEMC7
    /*
    AliAnalysisTaskBeautyCal *hfecalqaL07 = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqaL07);
    hfecalqaL07->SelectCollisionCandidates(AliVEvent::kEMC7);
    hfecalqaL07->SetElecIDsparse(FillElecSparse);
    hfecalqaL07->SetElecSys(FillElec);
    hfecalqaL07->SetTenderSwitch(UseTender);
    hfecalqaL07->SetThresholdEG1(thEG1ADC);
    hfecalqaL07->SetThresholdEG2(thEG2ADC);
    hfecalqaL07->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqaL07->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqaL07->SetCentralityMim(MimCent);
    hfecalqaL07->SetCentralityMax(MaxCent);
    hfecalqaL07->SetInvMassCut0(mass);
    hfecalqaL07->SetInvMassCut1(ptAss);
    hfecalqaL07->SetEtaRange(EtaRange);
    
    TString containerNameL07 = mgr->GetCommonFileName();
    containerNameL07 += ":PWGHF_hfeBeautyCalEMC7";
    containerNameL07 += ContNameExt;
    TString SubcontainerNameL07 = Form("BeautyCalEMC7_%s",calib);
    SubcontainerNameL07 += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerNameL07, TList::Class(),AliAnalysisManager::kOutputContainer, containerNameL07.Data());
    mgr->ConnectInput(hfecalqaL07, 0, cinput);
    mgr->ConnectOutput(hfecalqaL07, 1, coutput1);
    */
    /*
    // + kEMC8
    AliAnalysisTaskBeautyCal *hfecalqaL08 = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqaL08);
    hfecalqaL08->SelectCollisionCandidates(AliVEvent::kEMC8);
    hfecalqaL08->SetElecIDsparse(FillElecSparse);
    hfecalqaL08->SetTenderSwitch(UseTender);
    hfecalqaL08->SetThresholdEG1(thEG1ADC);
    hfecalqaL08->SetThresholdEG2(thEG2ADC);
    hfecalqaL08->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqaL08->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqaL08->SetCentralityMim(MimCent);
    hfecalqaL08->SetCentralityMax(MaxCent);
    
    TString containerNameL08 = mgr->GetCommonFileName();
    containerNameL08 += ":PWGHF_hfeBeautyCalEMC8";
    containerNameL08 += ContNameExt;
    TString SubcontainerNameL08 = Form("BeautyCalEMC8_%s",calib);
    SubcontainerNameL08 += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerNameL08, TList::Class(),AliAnalysisManager::kOutputContainer, containerNameL08.Data());
    mgr->ConnectInput(hfecalqaL08, 0, cinput);
    mgr->ConnectOutput(hfecalqaL08, 1, coutput1);
    
    */

    // + kcentral
    /*
    AliAnalysisTaskBeautyCal *hfecalqaCent = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqaCent);
    hfecalqaCent->SelectCollisionCandidates(AliVEvent::kCentral);
    hfecalqaCent->SetElecIDsparse(FillElecSparse);
    hfecalqaCent->SetElecSys(FillElec);
    hfecalqaCent->SetTenderSwitch(UseTender);
    hfecalqaCent->SetThresholdEG1(thEG1ADC);
    hfecalqaCent->SetThresholdEG2(thEG2ADC);
    hfecalqaCent->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqaCent->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqaCent->SetCentralityMim(MimCent);
    hfecalqaCent->SetCentralityMax(MaxCent);
    hfecalqaCent->SetInvMassCut0(mass);
    hfecalqaCent->SetInvMassCut1(ptAss);
    hfecalqaCent->SetEtaRange(EtaRange);
    
    TString containerNameCent = mgr->GetCommonFileName();
    containerNameCent += ":PWGHF_hfeBeautyCalCent";
    containerNameCent += ContNameExt;
    TString SubcontainerNameCent = Form("BeautyCalCent_%s",calib);
    SubcontainerNameCent += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerNameCent, TList::Class(),AliAnalysisManager::kOutputContainer, containerNameCent.Data());
    mgr->ConnectInput(hfecalqaCent, 0, cinput);
    mgr->ConnectOutput(hfecalqaCent, 1, coutput1);
    */
    
    // EMCal EGA
    if(ispPb)
    {
        // EMCal EGA EG1
        AliAnalysisTaskBeautyCal *hfecalqaTrig01 = new AliAnalysisTaskBeautyCal("emcqa");
        mgr->AddTask(hfecalqaTrig01);
        hfecalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig01->SetEMCalTriggerEG1(kTRUE);
        hfecalqaTrig01->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig01->SetElecSys(FillElec);
        hfecalqaTrig01->SetTenderSwitch(UseTender);
        hfecalqaTrig01->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig01->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig01->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaTrig01->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaTrig01->SetCentralityMim(MimCent);
        hfecalqaTrig01->SetCentralityMax(MaxCent);
        hfecalqaTrig01->SetInvMassCut0(mass);
        hfecalqaTrig01->SetInvMassCut1(ptAss);
        hfecalqaTrig01->SetEtaRange(EtaRange);
        hfecalqaTrig01->SetITSchi2(ITSchi2);
        hfecalqaTrig01->SetMinSig(nSig);
        hfecalqaTrig01->SetPileUpCut(pileup);
        
        TString containerName01 = mgr->GetCommonFileName();
        containerName01 += ":PWGHF_hfeBeautyCalTrigGAEG1";
        containerName01 += ContNameExt;
        TString SubcontainerName01 = Form("BeautyCalTrigGAEG1_%s",calib);
        SubcontainerName01 += ContNameExt;
        AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName01, TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
        mgr->ConnectInput(hfecalqaTrig01, 0, cinput);
        mgr->ConnectOutput(hfecalqaTrig01, 1, coutput1);
        
        // EMCal EGA EG2
        AliAnalysisTaskBeautyCal *hfecalqaTrig02 = new AliAnalysisTaskBeautyCal("emcqa");
        mgr->AddTask(hfecalqaTrig02);
        hfecalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig02->SetEMCalTriggerEG2(kTRUE);
        hfecalqaTrig02->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig02->SetElecSys(FillElec);
        hfecalqaTrig02->SetTenderSwitch(UseTender);
        hfecalqaTrig02->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig02->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig02->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaTrig02->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaTrig02->SetCentralityMim(MimCent);
        hfecalqaTrig02->SetCentralityMax(MaxCent);
        hfecalqaTrig02->SetInvMassCut0(mass);
        hfecalqaTrig02->SetInvMassCut1(ptAss);
        hfecalqaTrig02->SetEtaRange(EtaRange);
        hfecalqaTrig02->SetITSchi2(ITSchi2);        
        hfecalqaTrig02->SetMinSig(nSig);
        hfecalqaTrig02->SetPileUpCut(pileup);

        TString containerName02 = mgr->GetCommonFileName();
        containerName02 += ":PWGHF_hfeBeautyCalTrigGAEG2";
        containerName02 += ContNameExt;
        TString SubcontainerName02 = Form("BeautyCalTrigGAEG2_%s",calib);
        SubcontainerName02 += ContNameExt;
        AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName02, TList::Class(),AliAnalysisManager::kOutputContainer, containerName02.Data());
        mgr->ConnectInput(hfecalqaTrig02, 0, cinput);
        mgr->ConnectOutput(hfecalqaTrig02, 1, coutput1);
    }
    if(!ispPb)
    {
        // EMCal EGA
        AliAnalysisTaskBeautyCal *hfecalqaTrig0 = new AliAnalysisTaskBeautyCal("emcqa");
        mgr->AddTask(hfecalqaTrig0);
        hfecalqaTrig0->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        hfecalqaTrig0->SetElecIDsparse(FillElecSparse);
        hfecalqaTrig0->SetElecSys(FillElec);
        hfecalqaTrig0->SetTenderSwitch(UseTender);
        hfecalqaTrig0->SetThresholdEG1(thEG1ADC);
        hfecalqaTrig0->SetThresholdEG2(thEG2ADC);
        hfecalqaTrig0->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaTrig0->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaTrig0->SetCentralityMim(MimCent);
        hfecalqaTrig0->SetCentralityMax(MaxCent);
        hfecalqaTrig0->SetInvMassCut0(mass);
        hfecalqaTrig0->SetInvMassCut1(ptAss);
        hfecalqaTrig0->SetEtaRange(EtaRange);
        hfecalqaTrig0->SetITSchi2(ITSchi2);        
        hfecalqaTrig0->SetMinSig(nSig);
        hfecalqaTrig0->SetPileUpCut(pileup);

        TString containerName1 = mgr->GetCommonFileName();
        containerName1 += ":PWGHF_hfeBeautyCalTrigGA";
        containerName1 += ContNameExt;
        TString SubcontainerName1 = Form("BeautyCalTrigGAEG_%s",calib);
        SubcontainerName1 += ContNameExt;
        AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName1, TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
        mgr->ConnectInput(hfecalqaTrig0, 0, cinput);
        mgr->ConnectOutput(hfecalqaTrig0, 1, coutput1);
    }
    
    // MB
    AliAnalysisTaskBeautyCal *hfecalqaMB = new AliAnalysisTaskBeautyCal("emcqa");
    mgr->AddTask(hfecalqaMB);
    hfecalqaMB->SelectCollisionCandidates(AliVEvent::kMB);
    hfecalqaMB->SetElecIDsparse(FillElecSparse);
    hfecalqaMB->SetElecSys(FillElec);
    hfecalqaMB->SetTenderSwitch(UseTender);
    hfecalqaMB->SetThresholdEG1(thEG1ADC);
    hfecalqaMB->SetThresholdEG2(thEG2ADC);
    hfecalqaMB->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqaMB->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqaMB->SetCentralityMim(MimCent);
    hfecalqaMB->SetCentralityMax(MaxCent);
    hfecalqaMB->SetInvMassCut0(mass);
    hfecalqaMB->SetInvMassCut1(ptAss);
    hfecalqaMB->SetEtaRange(EtaRange);
    hfecalqaMB->SetITSchi2(ITSchi2);     
    hfecalqaMB->SetMinSig(nSig);
    hfecalqaMB->SetPileUpCut(pileup);     

    TString containerName2 = mgr->GetCommonFileName();
    containerName2 += ":PWGHF_hfeBeautyCalMB";
    containerName2 += ContNameExt;
    TString SubcontainerName2 = Form("BeautyCalMB_%s",calib);
    SubcontainerName2 += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName2, TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
    mgr->ConnectInput(hfecalqaMB, 0, cinput);
    mgr->ConnectOutput(hfecalqaMB, 1, coutput1); 
    
    //return hfecalqa;
    return NULL;
}
