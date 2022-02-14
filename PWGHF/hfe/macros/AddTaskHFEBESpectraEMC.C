AliAnalysisTask *AddTaskHFEBESpectraEMC(
                                 TString ContNameExt = "",
                                 Bool_t UseTender=kTRUE, TString CorrTaskClusCont = "caloClusters", TString CorrTaskTrkCont = "tracks",
                                 Bool_t FillElecSparse=kFALSE,
                                 Bool_t ClsTypeEMC=kTRUE, Bool_t ClsTypeDCAL=kTRUE,
                                 Bool_t SwitchPi0EtaWeightCalc = kTRUE,
                                 Bool_t SwitchNHFEeffi = kTRUE,
                                 Bool_t SwitchEleRecoEffi = kTRUE,
                                 Bool_t SwitchMCTempWeightCalc = kTRUE,
                                 Bool_t SwitchFillMCTemp = kTRUE,
                                 Bool_t SwitchRecalIP = kTRUE,
                                 Bool_t IsMC = kFALSE,
                                 Double_t deltaEta=0.05, Double_t deltaPhi=0.05,
                                 Double_t m02Min=0.05, Double_t m02Max1=0.9, Double_t m02Max2=0.9,
                                 Double_t m20Min=0.0, Double_t m20Max=20000,
                                 Double_t eovpMin=0.8, Double_t eovpMax=1.2,
                                 Int_t itsNCls = 3,
                                 Bool_t IsPPAnalysis=kFALSE,
                                 Bool_t SwitchEMCTrig = kFALSE,
                                 Bool_t hasTwoEMCTrigThres=kFALSE,
                                 Int_t MimCent = -1, Int_t MaxCent = -1,
                                 Bool_t IsPbPb2018 = kFALSE,
                                 Int_t thEG1ADC=140, Int_t thEG2ADC=89,
                                 TString centrality="V0M")
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

    //Add tender 
    char calib[100];
    if(UseTender)
    {
        if(MimCent==-1)
        {
            sprintf(calib,"wTend_");
        }
        else
        {
            sprintf(calib,"wTend_%d_%d",MimCent,MaxCent);
        }
        ContNameExt+=CorrTaskClusCont.Data();
    }
    else
    {
        if(MimCent==-1)
        {
            sprintf(calib,"woTend_");
        }
        else
        {
            sprintf(calib,"woTend_%d_%d",MimCent,MaxCent);
        }
    }
        
    if(ClsTypeEMC && !ClsTypeDCAL)ContNameExt+="_EMC";
    if(!ClsTypeEMC && ClsTypeDCAL)ContNameExt+="_DCAL";
    
    // INT7
    AliAnalysisTaskHFEBESpectraEMC *hfecalqa7 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
    mgr->AddTask(hfecalqa7);
    hfecalqa7->SelectCollisionCandidates(AliVEvent::kINT7);
    hfecalqa7->IsAnalysispp(IsPPAnalysis);
    hfecalqa7->IsMC(IsMC);
    hfecalqa7->SetElecIDsparse(FillElecSparse);
    hfecalqa7->SetTenderSwitch(UseTender);
    if(UseTender) hfecalqa7->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
    hfecalqa7->SetClusterTypeEMC(ClsTypeEMC);
    hfecalqa7->SetClusterTypeDCAL(ClsTypeDCAL);
    hfecalqa7->SetCentralityMim(MimCent);
    hfecalqa7->SetCentralityMax(MaxCent);
    hfecalqa7->SetCentralityEstimator(centrality.Data());
    hfecalqa7->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
    hfecalqa7->SetNonHFEEffi(SwitchNHFEeffi);
    hfecalqa7->SetElecRecoEffi(SwitchEleRecoEffi);
    hfecalqa7->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
    hfecalqa7->SwitchFillMCTemplate(SwitchFillMCTemp);
    hfecalqa7->SwitchRecalImpPar(SwitchRecalIP);
    hfecalqa7->SetTrackMatchPar(deltaEta, deltaPhi);
    hfecalqa7->SetM02Cut(m02Min,m02Max1,m02Max2);
    hfecalqa7->SetM20Cut(m20Min,m20Max);
    hfecalqa7->SetEovPCut(eovpMin,eovpMax);
    hfecalqa7->SetITSNCls(itsNCls);
    
    if(SwitchFillMCTemp){
        TString DMesonWeightMaps, BMesonWeightMaps;
        
        if(IsPPAnalysis){
            DMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/DMesonpTWeight.root";
            BMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/BMesonpTWeight.root";
        
            printf("\n### reading file %s ...\n",DMesonWeightMaps.Data());
            printf("\n### reading file %s ...\n",BMesonWeightMaps.Data());
        
            TFile* f2 = TFile::Open(DMesonWeightMaps.Data());
            if(f2){
                TH1 *D1 = (TH1*)f2->Get("RatD0");
                TH1 *D2 = (TH1*)f2->Get("RatD0Up");
                TH1 *D3 = (TH1*)f2->Get("RatD0Down");
            
                hfecalqa7->SetDmesonWeightHist(D1,D2,D3);
            }
            //  f2->Close();
            TFile* f3 = TFile::Open(BMesonWeightMaps.Data());
            if(f3){
                TH1 *B1 = (TH1*)f3->Get("RatBMes");
                TH1 *B2 = (TH1*)f3->Get("RatBMesMin");
                TH1 *B3 = (TH1*)f3->Get("RatBMesMax");
            
                hfecalqa7->SetBmesonWeightHist(B1,B2,B3);
            }
            //  f3->Close();
        }
        
        if(!IsPPAnalysis){
            DMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/CharmpTWeight_PbPb3050.root";
            BMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/BeautypTWeight_PbPb3050.root";
            
            printf("\n### reading file %s ...\n",DMesonWeightMaps.Data());
            printf("\n### reading file %s ...\n",BMesonWeightMaps.Data());
            
            TFile* f2 = TFile::Open(DMesonWeightMaps.Data());
            if(f2){
                TH1 *D0 = (TH1*)f2->Get("WeightD0");
                TH1 *DPlus = (TH1*)f2->Get("WeightDPlus");
                TH1 *Ds = (TH1*)f2->Get("WeightDs");
                TH1 *Lc = (TH1*)f2->Get("WeightLc");
                
                hfecalqa7->SetDmesonWeightHistPbPb(D0,DPlus,Ds,Lc);
            }
            TFile* f3 = TFile::Open(BMesonWeightMaps.Data());
            if(f3){
                TH1 *B = (TH1*)f3->Get("WeightB");
                
                hfecalqa7->SetBmesonWeightHistPbPb(B);
            }
        }
    }
    
    TString containerName7 = mgr->GetCommonFileName();
    containerName7 += ":PWGHF_HFEBESpectraEMC";
    containerName7 += ContNameExt;
    TString SubcontainerName7 = Form("HFEBESpectraEMC_INT7_%s",calib);
    SubcontainerName7 += ContNameExt;
    AliAnalysisDataContainer *cinput7  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(SubcontainerName7, TList::Class(),AliAnalysisManager::kOutputContainer, containerName7.Data());
    mgr->ConnectInput(hfecalqa7, 0, cinput7);
    mgr->ConnectOutput(hfecalqa7, 1, coutput7);
    
    if(IsPbPb2018){
        //Centrality trigger used in 2018
        
        AliAnalysisTaskHFEBESpectraEMC *hfecalqaCent = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
        mgr->AddTask(hfecalqaCent);
        if(MimCent == 0) hfecalqaCent->SelectCollisionCandidates(AliVEvent::kCentral);
        if(MimCent == 30) hfecalqaCent->SelectCollisionCandidates(AliVEvent::kSemiCentral);
        hfecalqaCent->IsAnalysispp(IsPPAnalysis);
        hfecalqaCent->IsMC(IsMC);
        hfecalqaCent->SetElecIDsparse(FillElecSparse);
        hfecalqaCent->SetTenderSwitch(UseTender);
        if(UseTender) hfecalqaCent->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
        hfecalqaCent->SetClusterTypeEMC(ClsTypeEMC);
        hfecalqaCent->SetClusterTypeDCAL(ClsTypeDCAL);
        hfecalqaCent->SetCentralityMim(MimCent);
        hfecalqaCent->SetCentralityMax(MaxCent);
        hfecalqaCent->SetCentralityEstimator(centrality.Data());
        hfecalqaCent->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
        hfecalqaCent->SetNonHFEEffi(SwitchNHFEeffi);
        hfecalqaCent->SetElecRecoEffi(SwitchEleRecoEffi);
        hfecalqaCent->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
        hfecalqaCent->SwitchFillMCTemplate(SwitchFillMCTemp);
        hfecalqaCent->SwitchRecalImpPar(SwitchRecalIP);
        hfecalqaCent->SetTrackMatchPar(deltaEta, deltaPhi);
        hfecalqaCent->SetM02Cut(m02Min,m02Max1,m02Max2);
        hfecalqaCent->SetM20Cut(m20Min,m20Max);
        hfecalqaCent->SetEovPCut(eovpMin,eovpMax);
        hfecalqaCent->SetITSNCls(itsNCls);
        
        TString containerNameCent = mgr->GetCommonFileName();
        containerNameCent += ":PWGHF_HFEBESpectraEMC_Cent";
        containerNameCent += ContNameExt;
        TString SubcontainerNameCent = Form("HFEBESpectraEMC_Cent_%s",calib);
        SubcontainerNameCent += ContNameExt;
        AliAnalysisDataContainer *cinputCent  = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutputCent = mgr->CreateContainer(SubcontainerNameCent, TList::Class(),AliAnalysisManager::kOutputContainer, containerNameCent.Data());
        mgr->ConnectInput(hfecalqaCent, 0, cinputCent);
        mgr->ConnectOutput(hfecalqaCent, 1, coutputCent);
    }
    
    if(SwitchEMCTrig)
    {
        // EMCal EGA EG1
        if(ClsTypeEMC){
            AliAnalysisTaskHFEBESpectraEMC *hfecalqaTrig01 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
            mgr->AddTask(hfecalqaTrig01);
            hfecalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
            hfecalqaTrig01->SetEMCalTriggerEG1(kTRUE);
            hfecalqaTrig01->IsAnalysispp(IsPPAnalysis);
            hfecalqaTrig01->IsMC(IsMC);
            hfecalqaTrig01->SetElecIDsparse(FillElecSparse);
            hfecalqaTrig01->SetTenderSwitch(UseTender);
            if(UseTender) hfecalqaTrig01->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
            hfecalqaTrig01->SetThresholdEG1(thEG1ADC);
            hfecalqaTrig01->SetThresholdEG2(thEG2ADC);
            hfecalqaTrig01->SetClusterTypeEMC(ClsTypeEMC);
            hfecalqaTrig01->SetClusterTypeDCAL(ClsTypeDCAL);
            hfecalqaTrig01->SetCentralityMim(MimCent);
            hfecalqaTrig01->SetCentralityMax(MaxCent);
            hfecalqaTrig01->SetCentralityEstimator(centrality.Data());
            hfecalqaTrig01->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
            hfecalqaTrig01->SetNonHFEEffi(SwitchNHFEeffi);
            hfecalqaTrig01->SetElecRecoEffi(SwitchEleRecoEffi);
            hfecalqaTrig01->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
            hfecalqaTrig01->SwitchFillMCTemplate(SwitchFillMCTemp);
            hfecalqaTrig01->SwitchRecalImpPar(SwitchRecalIP);
            hfecalqaTrig01->SetTrackMatchPar(deltaEta, deltaPhi);
            hfecalqaTrig01->SetM02Cut(m02Min,m02Max1,m02Max2);
            hfecalqaTrig01->SetM20Cut(m20Min,m20Max);
            hfecalqaTrig01->SetEovPCut(eovpMin,eovpMax);
            hfecalqaTrig01->SetITSNCls(itsNCls);
            
            TString containerName01 = mgr->GetCommonFileName();
            containerName01 += ":PWGHF_HFEBESpectraEMC_TrigGAEG1";
            containerName01 += ContNameExt;
            TString SubcontainerName01 = Form("HFEBESpectraEMC_EG1_%s",calib);
            SubcontainerName01 += ContNameExt;
            AliAnalysisDataContainer *cinput01  = mgr->GetCommonInputContainer();
            AliAnalysisDataContainer *coutput01 = mgr->CreateContainer(SubcontainerName01, TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
            mgr->ConnectInput(hfecalqaTrig01, 0, cinput01);
            mgr->ConnectOutput(hfecalqaTrig01, 1, coutput01);
        }
        
        if(ClsTypeDCAL){
            // DCal EGA DG1
            AliAnalysisTaskHFEBESpectraEMC *hfecalqaTrig01 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
            mgr->AddTask(hfecalqaTrig01);
            hfecalqaTrig01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
            hfecalqaTrig01->SetEMCalTriggerDG1(kTRUE);
            hfecalqaTrig01->IsAnalysispp(IsPPAnalysis);
            hfecalqaTrig01->IsMC(IsMC);
            hfecalqaTrig01->SetElecIDsparse(FillElecSparse);
            hfecalqaTrig01->SetTenderSwitch(UseTender);
            if(UseTender) hfecalqaTrig01->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
            hfecalqaTrig01->SetThresholdEG1(thEG1ADC);
            hfecalqaTrig01->SetThresholdEG2(thEG2ADC);
            hfecalqaTrig01->SetClusterTypeEMC(ClsTypeEMC);
            hfecalqaTrig01->SetClusterTypeDCAL(ClsTypeDCAL);
            hfecalqaTrig01->SetCentralityMim(MimCent);
            hfecalqaTrig01->SetCentralityMax(MaxCent);
            hfecalqaTrig01->SetCentralityEstimator(centrality.Data());
            hfecalqaTrig01->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
            hfecalqaTrig01->SetNonHFEEffi(SwitchNHFEeffi);
            hfecalqaTrig01->SetElecRecoEffi(SwitchEleRecoEffi);
            hfecalqaTrig01->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
            hfecalqaTrig01->SwitchFillMCTemplate(SwitchFillMCTemp);
            hfecalqaTrig01->SwitchRecalImpPar(SwitchRecalIP);
            hfecalqaTrig01->SetTrackMatchPar(deltaEta, deltaPhi);
            hfecalqaTrig01->SetM02Cut(m02Min,m02Max1,m02Max2);
            hfecalqaTrig01->SetM20Cut(m20Min,m20Max);
            hfecalqaTrig01->SetEovPCut(eovpMin,eovpMax);
            hfecalqaTrig01->SetITSNCls(itsNCls);
            
            TString containerName01 = mgr->GetCommonFileName();
            containerName01 += ":PWGHF_HFEBESpectraEMC_TrigGADG1";
            containerName01 += ContNameExt;
            TString SubcontainerName01 = Form("HFEBESpectraEMC_DG1_%s",calib);
            SubcontainerName01 += ContNameExt;
            AliAnalysisDataContainer *cinput01  = mgr->GetCommonInputContainer();
            AliAnalysisDataContainer *coutput01 = mgr->CreateContainer(SubcontainerName01, TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
            mgr->ConnectInput(hfecalqaTrig01, 0, cinput01);
            mgr->ConnectOutput(hfecalqaTrig01, 1, coutput01);
        }
    }
    
    if(SwitchEMCTrig && hasTwoEMCTrigThres)
    {
        if(ClsTypeEMC){
            // EMCal EGA EG2
            AliAnalysisTaskHFEBESpectraEMC *hfecalqaTrig02 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
            mgr->AddTask(hfecalqaTrig02);
            hfecalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
            hfecalqaTrig02->SetEMCalTriggerEG2(kTRUE);
            hfecalqaTrig02->IsAnalysispp(IsPPAnalysis);
            hfecalqaTrig02->IsMC(IsMC);
            hfecalqaTrig02->SetElecIDsparse(FillElecSparse);
            hfecalqaTrig02->SetTenderSwitch(UseTender);
            if(UseTender) hfecalqaTrig02->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
            hfecalqaTrig02->SetThresholdEG1(thEG1ADC);
            hfecalqaTrig02->SetThresholdEG2(thEG2ADC);
            hfecalqaTrig02->SetClusterTypeEMC(ClsTypeEMC);
            hfecalqaTrig02->SetClusterTypeDCAL(ClsTypeDCAL);
            hfecalqaTrig02->SetCentralityMim(MimCent);
            hfecalqaTrig02->SetCentralityMax(MaxCent);
            hfecalqaTrig02->SetCentralityEstimator(centrality.Data());
            hfecalqaTrig02->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
            hfecalqaTrig02->SetNonHFEEffi(SwitchNHFEeffi);
            hfecalqaTrig02->SetElecRecoEffi(SwitchEleRecoEffi);
            hfecalqaTrig02->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
            hfecalqaTrig02->SwitchFillMCTemplate(SwitchFillMCTemp);
            hfecalqaTrig02->SwitchRecalImpPar(SwitchRecalIP);
            hfecalqaTrig02->SetTrackMatchPar(deltaEta, deltaPhi);
            hfecalqaTrig02->SetM02Cut(m02Min,m02Max1,m02Max2);
            hfecalqaTrig02->SetM20Cut(m20Min,m20Max);
            hfecalqaTrig02->SetEovPCut(eovpMin,eovpMax);
            hfecalqaTrig02->SetITSNCls(itsNCls);

            TString containerName02 = mgr->GetCommonFileName();
            containerName02 += ":PWGHF_HFEBESpectraEMC_TrigGAEG2";
            containerName02 += ContNameExt;
            TString SubcontainerName02 = Form("HFEBESpectraEMC_EG2_%s",calib);
            SubcontainerName02 += ContNameExt;
            AliAnalysisDataContainer *cinput02  = mgr->GetCommonInputContainer();
            AliAnalysisDataContainer *coutput02 = mgr->CreateContainer(SubcontainerName02, TList::Class(),AliAnalysisManager::kOutputContainer, containerName02.Data());
            mgr->ConnectInput(hfecalqaTrig02, 0, cinput02);
            mgr->ConnectOutput(hfecalqaTrig02, 1, coutput02);
        }
        
        if(ClsTypeDCAL){
            // DCal EGA DG2
            AliAnalysisTaskHFEBESpectraEMC *hfecalqaTrig02 = new AliAnalysisTaskHFEBESpectraEMC("emcqa");
            mgr->AddTask(hfecalqaTrig02);
            hfecalqaTrig02->SelectCollisionCandidates(AliVEvent::kEMCEGA);
            hfecalqaTrig02->SetEMCalTriggerDG2(kTRUE);
            hfecalqaTrig02->IsAnalysispp(IsPPAnalysis);
            hfecalqaTrig02->IsMC(IsMC);
            hfecalqaTrig02->SetElecIDsparse(FillElecSparse);
            hfecalqaTrig02->SetTenderSwitch(UseTender);
            if(UseTender) hfecalqaTrig02->SetCorrectionTaskCont(CorrTaskClusCont, CorrTaskTrkCont);
            hfecalqaTrig02->SetThresholdEG1(thEG1ADC);
            hfecalqaTrig02->SetThresholdEG2(thEG2ADC);
            hfecalqaTrig02->SetClusterTypeEMC(ClsTypeEMC);
            hfecalqaTrig02->SetClusterTypeDCAL(ClsTypeDCAL);
            hfecalqaTrig02->SetCentralityMim(MimCent);
            hfecalqaTrig02->SetCentralityMax(MaxCent);
            hfecalqaTrig02->SetCentralityEstimator(centrality.Data());
            hfecalqaTrig02->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
            hfecalqaTrig02->SetNonHFEEffi(SwitchNHFEeffi);
            hfecalqaTrig02->SetElecRecoEffi(SwitchEleRecoEffi);
            hfecalqaTrig02->SwitchMCTemplateWeightCalc(SwitchMCTempWeightCalc);
            hfecalqaTrig02->SwitchFillMCTemplate(SwitchFillMCTemp);
            hfecalqaTrig02->SwitchRecalImpPar(SwitchRecalIP);
            hfecalqaTrig02->SetTrackMatchPar(deltaEta, deltaPhi);
            hfecalqaTrig02->SetM02Cut(m02Min,m02Max1,m02Max2);
            hfecalqaTrig02->SetM20Cut(m20Min,m20Max);
            hfecalqaTrig02->SetEovPCut(eovpMin,eovpMax);
            hfecalqaTrig02->SetITSNCls(itsNCls);
            
            TString containerName02 = mgr->GetCommonFileName();
            containerName02 += ":PWGHF_HFEBESpectraEMC_TrigGADG2";
            containerName02 += ContNameExt;
            TString SubcontainerName02 = Form("HFEBESpectraEMC_DG2_%s",calib);
            SubcontainerName02 += ContNameExt;
            AliAnalysisDataContainer *cinput02  = mgr->GetCommonInputContainer();
            AliAnalysisDataContainer *coutput02 = mgr->CreateContainer(SubcontainerName02, TList::Class(),AliAnalysisManager::kOutputContainer, containerName02.Data());
            mgr->ConnectInput(hfecalqaTrig02, 0, cinput02);
            mgr->ConnectOutput(hfecalqaTrig02, 1, coutput02);
        }
    }

    return NULL;
}

