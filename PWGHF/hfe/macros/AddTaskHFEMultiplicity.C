AliAnalysisTask* AddTaskHFEMultiplicity(TString suffixName = "",
                                        Bool_t readMC    =    kFALSE,
                                        Bool_t SwitchPi0EtaWeightCalc = kTRUE,
                                        Bool_t PhysSelINT7 =  kTRUE,
                                        Bool_t useTender   =  kTRUE,
                                        Bool_t ClsTypeEMC  =  kTRUE,
                                        Bool_t ClsTypeDCAL =  kTRUE,
                                        Bool_t flagEG1=  kTRUE ,
                                        Bool_t flagEG2=  kTRUE,
                                        Bool_t flagDG1=  kTRUE,
                                        Bool_t flagDG2=  kTRUE,
                                        TString estimatorFilename="",
                                        Bool_t zvtxcut =  kTRUE,
                                        Bool_t zvtxQA =  kTRUE,
                                        Double_t refMult=61.2,
                                        Int_t NcontV =2,
                                        Double_t RatioCrossedRowOverFindable=0.8,
                                        Int_t NTPCCrossRows = 100.,
                                        Int_t TPCNclus =80.,
                                        Int_t ITSNclus =3.,
                                        Double_t DeltaEta =0.01,
                                        Double_t DeltaPhi =0.01,
                                        Double_t DCAxyCut =2.4,
                                        Double_t DCAzCut =3.2,
                                        Double_t Etarange = 0.6,
                                        Double_t TrackPtMin = 2.,
                                        Double_t EopEMin = 0.8,
                                        Double_t EopEMax = 1.2,
                                        Double_t TPCNSigMin = -1.,
                                        Double_t TPCNSigMax = 3.,
                                        Double_t M02Min =0.02,
                                        Double_t M02Max =0.35,
                                        Int_t AssoTPCCluster = 80.,
                                        Int_t AssoITSCluster =3.,
                                        Double_t AssoEPt =0.1,
                                        Double_t AssoEEtarange =0.9,
                                        Double_t AssoENsigma =3.,
                                        Bool_t AssoITSRefit =kTRUE,
                                        Double_t InvmassCut = 0.14                    )

{ // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    Bool_t  fReadMC = kTRUE;
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
    if(!mcH){
        fReadMC=kFALSE;
    }
    
    
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";
    
    if(PhysSelINT7){
        AliAnalysisTaskHFEMultiplicity* HFEtaskINT7 = new AliAnalysisTaskHFEMultiplicity("");
        mgr->AddTask(HFEtaskINT7); //HFEtask is my task
        HFEtaskINT7->SetReadMC(readMC);
        HFEtaskINT7->SelectCollisionCandidates(AliVEvent::kINT7);
        HFEtaskINT7->SetTenderSwitch(useTender);
        HFEtaskINT7->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
        HFEtaskINT7->SetClusterTypeEMC(ClsTypeEMC);
        HFEtaskINT7->SetClusterTypeDCAL(ClsTypeDCAL);
        HFEtaskINT7->Setzvtxcut(zvtxcut);
        HFEtaskINT7->SetzvtxQA(zvtxQA);
        HFEtaskINT7->SetNcontVCut(NcontV);
        HFEtaskINT7->SetEtaRange(Etarange);
        HFEtaskINT7->SetRatioCrossedRowOverFindable(RatioCrossedRowOverFindable);
        HFEtaskINT7->SetNTPCCrossRows(NTPCCrossRows);
        HFEtaskINT7->SetNTPCCluster(TPCNclus);
        HFEtaskINT7->SetNITSCluster(ITSNclus);
        HFEtaskINT7->SetTrackpTMin(TrackPtMin);
        HFEtaskINT7->SetEMCalMatching(DeltaEta, DeltaPhi);
        HFEtaskINT7->SetDCACut(DCAxyCut,DCAzCut);
        HFEtaskINT7->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
        HFEtaskINT7->SetEopE(EopEMin,EopEMax);
        HFEtaskINT7->SetShowerShapeEM02(M02Min,M02Max);
        HFEtaskINT7->SetInvMassCut(InvmassCut);
        HFEtaskINT7->SetAssoTPCclus(AssoTPCCluster);
        HFEtaskINT7->SetAssoITSclus(AssoITSCluster);
        HFEtaskINT7->SetAssoITSrefit(AssoITSRefit);
        HFEtaskINT7->SetAssopTMin(AssoEPt);
        HFEtaskINT7->SetAssoEtarange(AssoEEtarange);
        HFEtaskINT7->SetAssoTPCnsig(AssoENsigma);
        
        TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
        if(!fileEstimator)  {
            Printf("File with multiplicity estimator not found\n");
            return NULL;
        }
        HFEtaskINT7->SetReferenceMultiplicity(refMult);
        const Char_t* profilebasename="SPDmult";
        
        const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
        TProfile* multEstimatorAvg[2];
        for(Int_t ip=0; ip<2; ip++) {
            cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
            multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
            if (!multEstimatorAvg[ip]) {
                Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
                return NULL;
            }
        }
        
        HFEtaskINT7->SetMultiProfileLHC16s(multEstimatorAvg[0]);
        HFEtaskINT7->SetMultiProfileLHC16r(multEstimatorAvg[1]);
        
        // Create containers for input/output
        TString finDirname         = "_INT7";
        TString outBasicname       = "EID";
        TString profname       = "coutputProf";
        
        finDirname           +=   suffixName.Data();
        outBasicname      +=   finDirname.Data();
        profname          +=   finDirname.Data();
        
        
        
        mgr->ConnectInput(HFEtaskINT7,0,mgr->GetCommonInputContainer());
        mgr->ConnectOutput(HFEtaskINT7,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        mgr->ConnectOutput(HFEtaskINT7,2,mgr->CreateContainer(profname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        
        
        
    }
    
    if(!PhysSelINT7){
        
        
        // EMCal EGA GA1
        
        AliAnalysisTaskHFEMultiplicity* HFEtaskGA = new AliAnalysisTaskHFEMultiplicity("");
        mgr->AddTask(HFEtaskGA);
        HFEtaskGA->SetReadMC(readMC);
        HFEtaskGA->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        HFEtaskGA->SetEMCalTriggerEG1(flagEG1);
        HFEtaskGA->SetEMCalTriggerDG1(flagDG1);
        HFEtaskGA->SetEMCalTriggerEG2(flagEG2);
        HFEtaskGA->SetEMCalTriggerDG2(flagDG2);
        HFEtaskGA->SetTenderSwitch(useTender);
        HFEtaskGA->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
        HFEtaskGA->SetClusterTypeEMC(ClsTypeEMC);
        HFEtaskGA->SetClusterTypeDCAL(ClsTypeDCAL);
        HFEtaskGA->Setzvtxcut(zvtxcut);
        HFEtaskGA->SetzvtxQA(zvtxQA);
        HFEtaskGA->SetNcontVCut(NcontV);
        HFEtaskGA->SetEtaRange(Etarange);
        HFEtaskGA->SetRatioCrossedRowOverFindable(RatioCrossedRowOverFindable);
        HFEtaskGA->SetNTPCCrossRows(NTPCCrossRows);
        HFEtaskGA->SetNTPCCluster(TPCNclus);
        HFEtaskGA->SetNITSCluster(ITSNclus);
        HFEtaskGA->SetTrackpTMin(TrackPtMin);
        HFEtaskGA->SetEMCalMatching(DeltaEta, DeltaPhi);
        HFEtaskGA->SetDCACut(DCAxyCut,DCAzCut);
        HFEtaskGA->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
        HFEtaskGA->SetEopE(EopEMin,EopEMax);
        HFEtaskGA->SetShowerShapeEM02(M02Min,M02Max);
        HFEtaskGA->SetInvMassCut(InvmassCut);
        HFEtaskGA->SetAssoTPCclus(AssoTPCCluster);
        HFEtaskGA->SetAssoITSclus(AssoITSCluster);
        HFEtaskGA->SetAssoITSrefit(AssoITSRefit);
        HFEtaskGA->SetAssopTMin(AssoEPt);
        HFEtaskGA->SetAssoEtarange(AssoEEtarange);
        HFEtaskGA->SetAssoTPCnsig(AssoENsigma);
        
        if(flagEG1 || flagDG1){
        
        TFile* fileEstimatorGA1=TFile::Open(estimatorFilename.Data());
        if(!fileEstimatorGA1)  {
            Printf("File with multiplicity estimator not found\n");
            return NULL;
        }
        HFEtaskGA->SetReferenceMultiplicity(refMult);
        const Char_t* profilebasenameEG1="SPDmultEG1";
        const Char_t* periodNamesGA1[2] = {"LHC16s", "LHC16r"};
        TProfile* multEstimatorAvgGA1[2];
        for(Int_t ip=0; ip<2; ip++) {
            cout<< " Trying to get "<<Form("%s_%s",profilebasenameEG1,periodNamesGA1[ip])<<endl;
            multEstimatorAvgGA1[ip] = (TProfile*)(fileEstimatorGA1->Get(Form("%s_%s",profilebasenameEG1,periodNamesGA1[ip]))->Clone(Form("%s_%s_clone",profilebasenameEG1,periodNamesGA1[ip])));
            if (!multEstimatorAvgGA1[ip]) {
                Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNamesGA1[ip]);
                return NULL;
            }
        }
        
        HFEtaskGA->SetMultiProfileLHC16s(multEstimatorAvgGA1[0]);
        HFEtaskGA->SetMultiProfileLHC16r(multEstimatorAvgGA1[1]);
        }
        if(flagEG2 || flagDG2){
        
            TFile* fileEstimatorGA2=TFile::Open(estimatorFilename.Data());
            if(!fileEstimatorGA2)  {
            Printf("File with multiplicity estimator not found\n");
            return NULL;
                    
            }
            HFEtaskGA->SetReferenceMultiplicity(refMult);
            const Char_t* profilebasenameEG2="SPDmultEG2";
            const Char_t* periodNamesGA2[2] = {"LHC16s", "LHC16r"};
            TProfile* multEstimatorAvgGA2[2];
            for(Int_t ip=0; ip<2; ip++) {
            cout<< " Trying to get "<<Form("%s_%s",profilebasenameEG2,periodNamesGA2[ip])<<endl;
            multEstimatorAvgGA2[ip] = (TProfile*)(fileEstimatorGA2->Get(Form("%s_%s",profilebasenameEG2,periodNamesGA2[ip]))->Clone(Form("%s_%s_clone",profilebasenameEG2,periodNamesGA2[ip])));
            if (!multEstimatorAvgGA2[ip]) {
                Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNamesGA2[ip]);
                return NULL;
                    }
                }
        
            HFEtaskGA->SetMultiProfileLHC16s(multEstimatorAvgGA2[0]);
            HFEtaskGA->SetMultiProfileLHC16r(multEstimatorAvgGA2[1]);
        
            }
        TString finDirnameGA         = "_GA";
        TString outBasicnameGA       = "EID";
        TString profnameGA       = "coutputProf";
        
        finDirnameGA           +=   suffixName.Data();
        outBasicnameGA      +=   finDirnameGA.Data();
        profnameGA          +=   finDirnameGA.Data();
        
        
        mgr->ConnectInput(HFEtaskGA,0,mgr->GetCommonInputContainer());
        mgr->ConnectOutput(HFEtaskGA,1,mgr->CreateContainer(outBasicnameGA, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        mgr->ConnectOutput(HFEtaskGA,2,mgr->CreateContainer(profnameGA, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        
    
    }
    
    
    
    return NULL;
}

