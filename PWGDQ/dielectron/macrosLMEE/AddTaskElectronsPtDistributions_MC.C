//____________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskElectronsPtDistributions_MC (
                                              
                                              Int_t       iset = 0,
                                              const char *file_CentralityBins,
                                              const char *file_Weights,
                                              Double_t    DCAxy_param0 = 0.00515869,
                                              Double_t    DCAxy_param1 = 0.01016680,
                                              Double_t    DCAxy_param2 = 1.34489,
                                              Double_t    DCAz_max     = 0.1,
                                              Int_t       ITS_minNcls     = 4,
                                              Int_t       TPC_minNcls     = 80,
                                              Int_t       TPC_nClsdEdx    = 50,
                                              Int_t       TPC_minCr       = 100,
                                              Double_t    MinCrOverFindableCls = 0.8,
                                              Double_t    MaxGoldenChi2 = 36,
                                              Double_t    MaxTPCchi2 = 4,
                                              Double_t    MaxITSchi2 = 3,
                                              Double_t    MaxFracSharedCls = 0.4,
                                              const char *ITSreq = "kFirst",
                                              Double_t    nsigmaTOF_max = 3.0,
                                              Double_t    nsigmaITS_max = 1.0,
                                              Double_t    nsigmaTPC_min = -3.0,
                                              Double_t    nsigmaTPC_max = 3.0,
                                              Double_t    massMin = 0.005,
                                              Double_t    massMax = 0.015,
                                              Double_t    phivLim = 40.0,
                                              Double_t    ptMin = 0.4,
                                              Double_t    ptMax = 5.0,
                                              Double_t    etaLim = 0.8
                                            )  {
   
    //Weights
    TFile *inputWeights = TFile::Open (file_Weights);
    TH2F *fHistoWeights[6];
    for ( Int_t ipart=0 ; ipart<6 ; ipart++ )
        fHistoWeights[ipart] = (TH2F*) inputWeights -> Get (Form("fHistoWeights(%d)",ipart));
    
    //Centrality Bins
    TFile *inputCentrBins = TFile::Open (file_CentralityBins);
    TH1F *fHistoCentrBins = (TH1F*) inputCentrBins->Get ("fHistoCentralityBins");

    
    //Load Libraries
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("ERROR: Analysis Manager Not Found!!!\n");
        return NULL;
    }
    
    
    //Retrieve Input Event Handler
    if (!mgr->GetInputEventHandler()) {
        printf("ERROR: No Input Event Handler!!!\n");
        return NULL;
    }
    
    
    //Input container
    AliAnalysisDataContainer *input = mgr->GetCommonInputContainer();

    
    //Analysis Task (Central)
    AliAnalysisTask *task_central = new AliAnalysisTask_Syst_PtDistributionsMC  (Form("PtDistributionsMC_Central(%d)",iset));
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetCentralityRange  (0.0,10.0);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetDCAparameters    (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetTrackCuts        (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,
                                                                                 MaxGoldenChi2,MaxTPCchi2,MaxITSchi2,MaxFracSharedCls,ITSreq);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetPIDCuts          (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetPrefilterCuts    (massMin,massMax,phivLim);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::SetKinematicCuts    (ptMin,ptMax,etaLim);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::GetCentralityBins           (fHistoCentrBins);
    task_central -> AliAnalysisTask_Syst_PtDistributionsMC::GetWeightsPtDistributions   (fHistoWeights[0], fHistoWeights[1],fHistoWeights[2],fHistoWeights[3],fHistoWeights[4],fHistoWeights[5]);
    mgr -> AddTask(task_central);
    AliAnalysisDataContainer *output_central = mgr->CreateContainer (Form("Input_Central_%d",iset),TList::Class(),AliAnalysisManager::kOutputContainer,"LMEEoutput.root");
    mgr -> ConnectInput  (task_central,0,input);
    mgr -> ConnectOutput (task_central,1,output_central);
    
    
    //Analysis Task (Semicentral)
    AliAnalysisTask *task_semicentral = new AliAnalysisTask_Syst_PtDistributionsMC  (Form("PtDistributionsMC_Semicentral(%d)",iset));
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetCentralityRange  (10.0,50.0);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetDCAparameters    (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetTrackCuts        (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,
                                                                                     MaxGoldenChi2,MaxTPCchi2,MaxITSchi2,MaxFracSharedCls,ITSreq);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetPIDCuts          (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetPrefilterCuts    (massMin,massMax,phivLim);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::SetKinematicCuts    (ptMin,ptMax,etaLim);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::GetCentralityBins           (fHistoCentrBins);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsMC::GetWeightsPtDistributions   (fHistoWeights[0], fHistoWeights[1],fHistoWeights[2],fHistoWeights[3],fHistoWeights[4],fHistoWeights[5]);
    mgr -> AddTask(task_semicentral);
    AliAnalysisDataContainer *output_semicentral = mgr->CreateContainer (Form("Input_Semicentral_%d",iset),TList::Class(),AliAnalysisManager::kOutputContainer,"LMEEoutput.root");
    mgr -> ConnectInput  (task_semicentral,0,input);
    mgr -> ConnectOutput (task_semicentral,1,output_semicentral);
    
    return task;
}
//____________________________________________________________________________________________________________________________________________________________________________________________________________

