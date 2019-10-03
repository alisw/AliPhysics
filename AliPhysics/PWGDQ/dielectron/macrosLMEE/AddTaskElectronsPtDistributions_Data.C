//____________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskElectronsPtDistributions_Data (
                                              
                                              Int_t       iset = 0,
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
    AliAnalysisTask *task_central = new AliAnalysisTask_Syst_PtDistributionsData  (Form("PtDistributionsData_Central(%d)",iset));
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetCentralityRange  (0.0,10.0);
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetDCAparameters    (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetTrackCuts        (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,
                                                                                   MaxGoldenChi2,MaxTPCchi2,MaxITSchi2,MaxFracSharedCls,ITSreq);
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetPIDCuts          (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetPrefilterCuts    (massMin,massMax,phivLim);
    task_central -> AliAnalysisTask_Syst_PtDistributionsData::SetKinematicCuts    (ptMin,ptMax,etaLim);
    mgr -> AddTask (task_central);
    AliAnalysisDataContainer *output_central = mgr->CreateContainer (Form("Input_Central_%d",iset),TList::Class(),AliAnalysisManager::kOutputContainer,"LMEEoutput.root");
    mgr -> ConnectInput  (task_central,0,input);
    mgr -> ConnectOutput (task_central,1,output_central);
    
    
    //Analysis Task (Semicentral)
    AliAnalysisTask *task_semicentral = new AliAnalysisTask_Syst_PtDistributionsData  (Form("PtDistributionsData_Semicentral(%d)",iset));
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetCentralityRange  (10.0,50.0);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetDCAparameters    (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetTrackCuts        (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,
                                                                                       MaxGoldenChi2,MaxTPCchi2,MaxITSchi2,MaxFracSharedCls,ITSreq);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetPIDCuts          (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetPrefilterCuts    (massMin,massMax,phivLim);
    task_semicentral -> AliAnalysisTask_Syst_PtDistributionsData::SetKinematicCuts    (ptMin,ptMax,etaLim);
    mgr -> AddTask(task_semicentral);
    AliAnalysisDataContainer *output_semicentral = mgr->CreateContainer (Form("Input_Semicentral_%d",iset),TList::Class(),AliAnalysisManager::kOutputContainer,"LMEEoutput.root");
    mgr -> ConnectInput  (task_semicentral,0,input);
    mgr -> ConnectOutput (task_semicentral,1,output_semicentral);
    
    
    return task;
}
//____________________________________________________________________________________________________________________________________________________________________________________________________________

