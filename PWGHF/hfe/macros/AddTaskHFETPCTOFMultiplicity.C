class AliAnalysisDataContainer;
AliAnalysisTask* AddTaskHFETPCTOFMultiplicity(TString suffixName = "",
                                              Bool_t readMC    =    kFALSE,
                                              Bool_t SwitchPi0EtaWeightCalc = kTRUE,
                                              Bool_t PhysSelINT7 =  kTRUE,
                                              TString estimatorFilename="",
                                              Double_t refMult=61.2,
                                              Int_t NcontV =2,
                                              Int_t MaxTPCclus = 100.,
                                              Int_t TPCNclus =80.,
                                              Int_t ITSNclus =3.,
                                              Bool_t SPDBoth= kTRUE ,
                                              Bool_t SPDAny= kFALSE ,
                                              Bool_t SPDFirst= kFALSE ,
                                              Double_t DCAxyCut =2.4,
                                              Double_t DCAzCut =3.2,
                                              Double_t Etarange = 0.7,
                                              Double_t CutNsigmaTOF = 3.,
                                              Double_t TPCNSigMin = -1.,
                                              Double_t TPCNSigMax = 3.,
                                              Int_t AssoTPCCluster = 80.,
                                              Int_t AssoITSCluster =3.,
                                              Double_t AssoEPt =0.1,
                                              Double_t AssoEEtarange =0.9,
                                              Double_t AssoENsigma =3.,
                                              Bool_t AssoITSRefit =kTRUE,
                                              Double_t InvmassCut = 0.14)

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
        AliAnalysisTaskHFETPCTOFMultiplicity* HFEtaskINT7 = new AliAnalysisTaskHFETPCTOFMultiplicity("");
        mgr->AddTask(HFEtaskINT7); //HFEtask is my task
        HFEtaskINT7->SetReadMC(readMC);
        HFEtaskINT7->SelectCollisionCandidates(AliVEvent::kINT7);
        HFEtaskINT7->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeightCalc);
        HFEtaskINT7->SetNcontVCut(NcontV);
        HFEtaskINT7->SetEtaRange(Etarange);
        HFEtaskINT7->SetMaxTPCCluster(MaxTPCclus);
        HFEtaskINT7->SetNTPCCluster(TPCNclus);
        HFEtaskINT7->SetNITSCluster(ITSNclus);
        HFEtaskINT7->SetTOFNSigmaCut(CutNsigmaTOF);
        HFEtaskINT7->SetDCACut(DCAxyCut,DCAzCut);
        HFEtaskINT7->SetTPCnsigma(TPCNSigMin,TPCNSigMax);
        HFEtaskINT7->SetInvMassCut(InvmassCut);
        HFEtaskINT7->SetAssoTPCclus(AssoTPCCluster);
        HFEtaskINT7->SetAssoITSclus(AssoITSCluster);
        HFEtaskINT7->SetAssoITSrefit(AssoITSRefit);
        HFEtaskINT7->SetAssopTMin(AssoEPt);
        HFEtaskINT7->SetAssoEtarange(AssoEEtarange);
        HFEtaskINT7->SetAssoTPCnsig(AssoENsigma);
        HFEtaskINT7->SetHitsOnSPDLayers(SPDBoth,SPDAny,SPDFirst);
        TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
        if(!fileEstimator)  {
            printf("File with multiplicity estimator not found\n");
            return 0x0;
        }
        HFEtaskINT7->SetReferenceMultiplicity(refMult);
        const Char_t* profilebasename="SPDmult";
        
        const Char_t* periodNames[2] = {"LHC16s", "LHC16r"};
        TProfile* multEstimatorAvg[2];
        for(Int_t ip=0; ip<2; ip++) {
            cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
            multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
            if (!multEstimatorAvg[ip]) {
                printf("Multiplicity estimator not found! Please check your estimator file");
                return 0x0;
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
    
    return NULL;
    
}

