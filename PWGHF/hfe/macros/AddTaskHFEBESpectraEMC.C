AliAnalysisTask *AddTaskHFEBESpectraEMC(
                                 Bool_t UseTender=kTRUE,
                                 Bool_t FillElecSparse=kFALSE,
                                 Bool_t ClsTypeEMC=kTRUE, Bool_t ClsTypeDCAL=kTRUE,
                                 Bool_t SwitchPi0EtaWeightCalc = kTRUE,
                                 Bool_t SwitchNHFEeffi = kTRUE,
                                 Int_t MimCent = -1, Int_t MaxCent = -1,
                                 TString ContNameExt = "", TString centrality="V0M",
                                 Bool_t hasTwoEMCTrigThres=kFALSE, Int_t thEG1ADC=140, Int_t thEG2ADC=89)
{
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTask", "No analysis manager found.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTask", "This task requires an input event handler");
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
    
    // INT7
    AliAnalysisTaskHFEBESpectraEMC *hfecalqa7 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
    mgr->AddTask(hfecalqa7);
    hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
    hfecalqa7->SetElecIDsparse(FillElecSparse);
    hfecalqa7->SetTenderSwitch(UseTender);
    hfecalqa7->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqa7->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqa7->SetCentralityMim(MimCent);
    hfecalqa7->SetCentralityMax(MaxCent);
    hfecalqa7->SetCentralityEstimator(centrality.Data());
    hfecalqa7->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
    hfecalqa7->SetNonHFEEffi(SwitchNHFEeffi);
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PWGHF_HFEBESpectraEMC";
    containerName7 += ContNameExt;
    TString SubcontainerName7 = Form("HFEBESpectraEMC_INT7_%s",calib);
    SubcontainerName7 += ContNameExt;
    AliAnalysisDataContainer *cinput7  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(SubcontainerName7, TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    mgr->ConnectInput(hfecalqa7, 0, cinput7);
    mgr->ConnectOutput(hfecalqa7, 1, coutput7);
    
    return NULL;
}

