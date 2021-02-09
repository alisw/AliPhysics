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
                                                         TString ElectronEfficiencyFile = "",
                                                         TString BackgroundWFileToData = "alien:///alice/cern.ch/user/h/hzanoli/BackgroundW/BackgroundWToData.root",
                                                         TString Sufix = "",
                                                         Float_t ZvtxMin  = -10.,
                                                         Float_t ZvtxMax = 10.,
                                                         Int_t ZvtxBinConfig = 0,
                                                         Bool_t UseTOF = kTRUE,
                                                         Int_t NBinsCentrality = 4
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
    taskName += Sufix;
    
    AliAnalysisTaskHFEpACorrelation *task = ConfigHFEpACorrelation(taskName, Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID,UseGlobalTracksForHadrons,CentralityEstimator,UseTOF);
    
    
    Double_t vertexBins0[] = {-10, -5.3, -2.6, -0.5, 1.5, 3.5, 5.7,10.01};
    Double_t vertexBins1[] = {-10, -5, -3, -1, 1, 3, 5, 10.01};
    Double_t vertexBins2[] = {-10, -6, -2.2, -1, 0.8, 2.7, 5.2,10.01};
    Double_t vertexBins3[] = {-10, -4, -1.0,0.7, 1.9, 3.7, 6.2,10.01};
    Double_t vertexBins4[] = {-10, -5.5, -3.0,-0.9, 1.2, 2.9, 5.9,10.01};
    
    switch (ZvtxBinConfig) {
        case 1:
            task->SetZvtxBins(8,vertexBins1);
            break;
        case 2:
            task->SetZvtxBins(8,vertexBins2);
            break;
        case 3:
            task->SetZvtxBins(8,vertexBins3);
            break;
        case 4:
            task->SetZvtxBins(8,vertexBins4);
            break;
            
        default:
            task->SetZvtxBins(8,vertexBins0);
            break;
    }
    
    TArrayD CentralityArray(NBinsCentrality+1);
    
    for (Int_t i = 0; i <= NBinsCentrality; i++)
    {
        CentralityArray.SetAt(MultiplicityLow + (Double_t)i*(MultiplicityUp-MultiplicityLow)/NBinsCentrality, i);
    }
    
    task->SetCentralityBins(NBinsCentrality+1,CentralityArray.GetArray());

    
    //_______________________
    //Trigger
    
    if (!ispp)
        task->SelectCollisionCandidates(AliVEvent::kINT7);
    else
        task->SelectCollisionCandidates(AliVEvent::kMB);
    
    if(pTBin ==0)
    {
        Float_t pTBinsCorrelation[] = {0.5,0.75,1.0,1.25,1.5,2,2.5,3,4,6};
        task->SetpTBins(10,pTBinsCorrelation);
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
    
    task->SetZVtxCut(ZvtxMin,ZvtxMax);

    
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
        
        TFile *ElectronEffFile_f =  TFile::Open(ElectronEfficiencyFile.Data());
        if (ElectronEffFile_f)
        {
            
            TH3F* test = (TH3F*) ElectronEffFile_f->Get("EffE");
            if (test)
                task->SetEfficiencyElectron(test);
        }

        
        
    }
    
    if (isMC)
    {
        
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
    
    
   
    
    TString containerName;
    containerName += "HFE_h";
    containerName += Form("_%d_%d_%d_%d_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_%d_%d_%d",pTBin,Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID);
    containerName += Sufix.Data();
    
    TString fileName = mgr->GetCommonFileName();
    fileName.Append(":");
    fileName += containerName;
    
    //Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(containerName.Data(), TList::Class(),    AliAnalysisManager::kOutputContainer, fileName.Data());
    
    //Connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);
    
    return task;
}
