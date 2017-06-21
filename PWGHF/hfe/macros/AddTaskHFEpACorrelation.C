AliAnalysisTaskHFEpACorrelation *AddTaskHFEpACorrelation(
                                                         Int_t pTBin = 0,
                                                         Bool_t Correlation = kTRUE,
                                                         Bool_t ispp = kFALSE,
                                                         Bool_t isMC = kTRUE,
                                                         Double_t ElectronDCAxy = 0.25,
                                                         Double_t ElectronDCAz = 1.0,
                                                         Double_t HadronDCAxy = 0.25,
                                                         Double_t HadronDCAz = 1.0,
                                                         Double_t TPCPIDLow = -0.5,
                                                         Double_t TPCPIDUp = 3.0,
                                                         Double_t InvariantMassCut = 0.14,
                                                         Double_t pTCutPartner = 0.0,
                                                         Double_t MultiplicityLow = 0.,
                                                         Double_t MultiplicityUp = 100.,
                                                         Double_t HadronPtCutLow = 0.3,
                                                         Double_t HadronPtCutUp = 2.0,
                                                         Double_t EtaCutLow = -0.8,
                                                         Double_t EtaCutUp = 0.8,
                                                         Double_t NonHFEangleCut = 999,
                                                         Int_t NHitsITS = 4,
                                                         Int_t SPDLayers = 0,
                                                         Int_t TPCNCluster = 100,
                                                         Int_t TPCNClusterPartner = 60,
                                                         Int_t TPCNClusterPID = 80,
                                                         Bool_t UseGlobalTracksForHadrons = kTRUE,
                                                         Int_t CentralityEstimator = 0,
                                                         TString HadronEfficiencyFile = "alien:///alice/cern.ch/user/h/hzanoli/Efficiency/Hadron_Tracking.root",
                                                         TString BackgroundWFile = "alien:///alice/cern.ch/user/h/hzanoli/BackgroundW/BackgroundW.root",
                                                         TString BackgroundWFileToData = "alien:///alice/cern.ch/user/h/hzanoli/BackgroundW/BackgroundWToData.root"
                                                         )
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    if (!mgr) {
        ::Error("AddTaskEMCalHFEpA", "No analysis manager to connect to.");
        return NULL;
    }
    
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskEMCalHFEpA", "This task requires an input event handler");
        return NULL;
    }
    
    //_______________________
    //Config Task
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEpACorrelation.C");
    TString taskName = "HFe_h";
    taskName.Append(Form("%d_%d_%d_%d_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_%d_%d_%d_%d_%d",pTBin, Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID,UseGlobalTracksForHadrons,CentralityEstimator));
    
    AliAnalysisTaskHFEpACorrelation *task = ConfigHFEpACorrelation(taskName, Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID,UseGlobalTracksForHadrons,CentralityEstimator);
    
    //_______________________
    //Trigger
    
    if (!ispp)
        task->SelectCollisionCandidates(AliVEvent::kINT7);
    else
        task->SelectCollisionCandidates(AliVEvent::kMB);
    
    if(pTBin ==0)
    {
        Float_t pTBinsCorrelation[] = {0.5,0.75,1.0,1.25,1.5,2,2.5,3,4,6};
        task->SetpTBins(11,pTBinsCorrelation);
    }
    else if (pTBin ==1)
    {
        Float_t pTBinsCorrelation[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4, 5, 6}
        task->SetpTBins(12,pTBinsCorrelation);
    }
    
    if(pTBin ==2)
    {
        Float_t pTBinsCorrelation[] = {0.5,1,2,4,6};
        task->SetpTBins(5,pTBinsCorrelation);
    }
    
    else if (pTBin ==3)
    {
        Float_t pTBinsCorrelation[] = {0.5, 0.7, 0.9, 1.0, 1.2, 1.4, 1.6,1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6}
        task->SetpTBins(16,pTBinsCorrelation);
    }
    
    
    if(Correlation)
    {
        TFile *fileH = TFile::Open(HadronEfficiencyFile.Data());
        TH3F* HadronEffHisto = (TH3F*) fileH->Get(Form("HadronEff_%1.2f_%1.2f",HadronDCAxy,HadronDCAz));
        if(HadronEffHisto)
            task->SetEfficiencyHadron(HadronEffHisto);
        else
        {
            printf("=!=!=!=!=!=!=!=! No Hadron Correction =!=!=!=!=!=!=!=!\n");
            printf("=!=!=!=!=!=!=!=! Correlation will be useless =!=!=!=!=!=!=!=!\n");
        }
    }
    
    if (isMC)
    {
        TFile *fileBkgW= TFile::Open(BackgroundWFile.Data());
        if (fileBkgW)
        {
            TH1F *Pi0 = (TH1F*) fileBkgW->Get("Pi0W");
            TH1F *Eta = (TH1F*) fileBkgW->Get("EtaW");
            
            if (Pi0)
                task->SetBackgroundPi0Weight(Pi0);
            else
                printf("MC analysis with no pi0 weight hijing\n");
            
            if (Eta)
                task->SetBackgroundEtaWeight(Eta);
            else
                printf("MC analysis with no Eta weight hijing \n");
            
        }
        else
            printf("Background weight to Hijing not available\n");
        
        TFile *fileBkgWToData= TFile::Open(BackgroundWFileToData.Data());
        
        if (fileBkgWToData)
        {
            TH1F *Pi0W = (TH1F*) fileBkgWToData->Get("Pi0");
            TH1F *EtaW = (TH1F*) fileBkgWToData->Get("Eta");
            
            if (Pi0W)
                task->SetBackgroundPi0WeightToData(Pi0W);
            else
                printf("MC analysis with no pi0 weight to data \n");
            
            if (Eta)
                task->SetBackgroundEtaWeightToData(EtaW);
            else
                printf("MC analysis with no Eta weight to data\n");
            
        }
        else
            printf("Background weight to Data not available. Check the path to the file! \n");

    }
    
    
    TString containerName = mgr->GetCommonFileName();
    containerName += ":HFE_h";
    containerName += Form("_%d_%d_%d_%d_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_%d_%d_%d",pTBin,Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID);
    
    //Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer( Form("eh_%d_%d_%d_%d_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_%d_%d_%d",pTBin,Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());
    
    //Connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);
    
    return task;
}
