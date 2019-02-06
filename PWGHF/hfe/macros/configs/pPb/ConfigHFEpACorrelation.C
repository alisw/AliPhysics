/* Possible Configurations
 
 Electron DCA cut
 ElectronDCAxy, ElectronDCAz
 
 ITS Number of Hits = NHitsITS (obvius)
 
 SPD
 0 = kBoth
 1 = kFirst
 2 = kAny
 
 Number of Clusters On TPC  = TPCNCluster (obvius)
 
 TPC PID Cluster = TPCNClusterPID (obvius)
 TPCNSigma = TPCNClusterPID (obvius)
 
CentralityEstimator 
0 = ZNA (Standard)
1 = VOA (used before)

 
 */

AliAnalysisTaskHFEpACorrelation* ConfigHFEpACorrelation(TString taskName = "HFe_h",
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
                                                        Bool_t UseTOF = kTRUE)
{
    
    
    
    
    
    ///_______________________________________________________________________________________________________________
    ///Track selection: Cuts used to ensure a minimum quality level of the tracks selected to perform the analysis
    
    AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsMinBias","HFE Cuts");
    hfecuts->CreateStandardCuts();
    
    printf("============ Configuring AliHFEcuts ============ \n");
    hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetMinNClustersTPC(TPCNCluster);
    hfecuts->SetMinNClustersTPCPID(TPCNClusterPID); 						                    //Minimum number of clusters for dE/dx
    printf("TPCNCluster = %d and TPCNClusterPID = %d \n",TPCNCluster,TPCNClusterPID);
    hfecuts->SetMinRatioTPCclusters(0.6);
    
    //ITS
    printf("ITS pixel = %d \n", SPDLayers);
    if (SPDLayers == 0)
    {
        hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
        printf("ITS pixel =  kBoth\n");
    }
    else if (SPDLayers == 1)
    {
        hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
        printf("ITS pixel =  kFirst\n");
    }
    else if (SPDLayers == 2)
    {
        hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
        printf("ITS pixel =  kAny\n");
    }
    else
    {
        hfecuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
        printf("ERROR: ITSPIXEL not set. Setting it as ITS pixel =  kBoth\n");
    }
    
    hfecuts->SetMinNClustersITS(NHitsITS);
    printf("Number of Points on the ITS: %d\n", NHitsITS);
    
    hfecuts->SetCheckITSLayerStatus(kFALSE);
    
    //Additional Cuts
    hfecuts->SetPtRange(0.5, 6);//Transversal momentum range in GeV/c
    
    hfecuts->SetMaxImpactParam(ElectronDCAxy,ElectronDCAz);
    //hfecuts->SetUseMixedVertex(kTRUE); //Copy Jan
    printf("Electron DCA cut to xy = %1.2fcm z = %1.2f \n", ElectronDCAxy,ElectronDCAz);
    
    
    //Event Selection
    hfecuts->SetVertexRange(10.);													//
    //hfecuts->SetProductionVertex(0,0.3,0,0.3);									//
    ///_______________________________________________________________________________________________________________
    
    ///_________________________________________________________________________________________________________________________
    ///Task config
    AliAnalysisTaskHFEpACorrelation *task = new AliAnalysisTaskHFEpACorrelation(taskName.Data());
    printf("task ------------------------ %p\n ", task);
    
    task->SetHFECuts(hfecuts);
    if(Correlation)
    {
        task->SetCorrelationAnalysis();
        task->SetEventMixing(kTRUE);
        printf("Correlation Analysis ON \n");
    }
    task->SetAODanalysis(kTRUE);
    if(ispp)
        task->SetPPanalysis(kTRUE);
    
    //On data it is defined on the HFEcuts, but we should use the same cuts for the NHF partner
    task->SetdcaCut(ElectronDCAxy,ElectronDCAz);
    if (UseGlobalTracksForHadrons)
        task->UseGlobalTracksHadron();
    task->SetAssHadronPtRange(HadronPtCutLow,HadronPtCutUp);
    printf("HadronPtCutLow = %1.2f HadronPtCutUp = %1.2f, Using GlobalTracks = %d\n", HadronPtCutLow,HadronPtCutUp,UseGlobalTracksForHadrons);
    
    
    task->SetAdditionalCuts(pTCutPartner,TPCNClusterPartner);
    printf("pTCutPartner = %1.2f TPCNClusterPartner = %1.2f\n", pTCutPartner,TPCNClusterPartner);
    
    task->SetNonHFEmassCut(InvariantMassCut);
    printf("Invariant mass cut = %1.2f", InvariantMassCut);
    
    task->SetEtaCut(EtaCutLow,EtaCutUp);
    printf("EtaCutLow = %1.2f EtaCutUp = %1.2f\n",EtaCutLow,EtaCutUp);
    
    task->SetNonHFEangleCut(NonHFEangleCut);
    printf("NonHFEangleCut = %1.2f\n",NonHFEangleCut);
    
    
    task->SetUseDCACutHadron();
    task->SetDCACutHadron(HadronDCAxy, HadronDCAz);
    printf("Hadron DCA cut to xy = %1.2fcm z = %1.2f \n", HadronDCAxy, HadronDCAz);
    
    
    task->SetCentrality(MultiplicityLow,MultiplicityUp);
    printf("MinMultiplicity = %1.2f MaxMultiplicy = %1.2f\n", MultiplicityLow,MultiplicityUp);
    
    task->SetCentralityEstimator(CentralityEstimator);
    printf("Centrality estimator = %d",CentralityEstimator);
    
    
    ///_______________________________________________________________________________________________________________
    
    ///_______________________________________________________________________________________________________________
    ///Particle identification
    AliHFEpid *pid = task->GetPID();
    
    //______________________________________
    //In the case of a simulation
    if(isMC)
    {
        pid->SetHasMCData(kTRUE);
        task->SetMCanalysis();
    }
    //______________________________________
    
    //______________________________________________________
    //Configure PID
    //_________________________
    //TPC+TOF PID
    
    if (UseTOF)
    {
        pid->AddDetector("TOF", 0);				//Add TOF PID
        pid->AddDetector("TPC", 1);				//Add TPC PID
    }
    else
        pid->AddDetector("TPC", 0);				//Add TPC PID

    
    //_________________________
    //Configure TPC cut
    //Defaul = -1 to 3 sigmas
    //Note that it is also possible to define a model instead of a constant
    //--------->For this change the "cut model"
    
    Double_t params[4];
    char *cutmodel;
    cutmodel = "pol0";
    
    params[0] = TPCPIDLow;
    
    pid->ConfigureTPCdefaultCut(cutmodel,params,TPCPIDUp);
    //_______________________________________________________
    ///_______________________________________________________________________________________________________________
    
    printf("*************************************\n");
    printf("Configuring Task:\n");
    pid->PrintStatus();
    printf("*************************************\n");
    
    return task;
}
