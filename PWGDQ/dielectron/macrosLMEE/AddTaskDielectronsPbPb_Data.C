//____________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskDielectronsPbPb_Data (
                                              
                                              Int_t    iset = 0,
                                              const char *centrality = "Central",
                                              Double_t CentralityMin =  0.0,
                                              Double_t CentralityMax = 10.0,
                                              const char *magfield = "up",
                                              Int_t    MaxNumberEvts = 10000,
                                              Int_t    MaxNumberTrks = 100000,
                                              Int_t    NumberEvtsToMix = 20,
                                              Int_t    NcentralityBins = 5,
                                              Int_t    NvertexBins = 10,
                                              Int_t    NeventPlaneBins = 1,
                                              Double_t DCAxy_param0 = 0.00515869,
                                              Double_t DCAxy_param1 = 0.01016680,
                                              Double_t DCAxy_param2 = 1.34489,
                                              Double_t DCAz_max     = 0.1,
                                              Int_t    ITS_minNcls     = 4,
                                              Int_t    TPC_minNcls     = 80,
                                              Int_t    TPC_nClsdEdx    = 50,
                                              Int_t    TPC_minCr       = 100,
                                              Double_t MinCrOverFindableCls = 0.8,
                                              Double_t MaxGoldenChi2 = 36,
                                              Double_t MaxTPCchi2 = 4,
                                              Double_t MaxITSchi2 = 3,
                                              Double_t MaxFracSharedCls = 0.4,
                                              const char *ITSreq = "kFirst",
                                              Double_t nsigmaTOF_max = 3.0,
                                              Double_t nsigmaITS_max = 1.0,
                                              Double_t nsigmaTPC_min = -3.0,
                                              Double_t nsigmaTPC_max = 3.0,
                                              Double_t massLim = 0.02,
                                              Double_t phivLim = 60.0,
                                              Double_t ptMin = 0.4,
                                              Double_t ptMax = 5.0,
                                              Double_t etaLim = 0.8
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
    
    
    //Task Name,InputBox & Output File
    const char *TaskName =   Form("DielectronTask_Data_Set%d_%s_MagField_%s",iset,centrality,magfield);
    const char *InputBox =   Form("Input_%d_%s_%s",iset,centrality,magfield);
    const char *OutputFile = Form("InputFile_DielectronsPbPb_Data_Set%d_%s_MagField_%s.root",iset,centrality,magfield);
    
    
    //Input container
    AliAnalysisDataContainer *input = mgr->GetCommonInputContainer();

    
    //Analysis Task
    AliAnalysisTask *task = new AliAnalysisTaskDielectronsPbPb_Data  (TaskName);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetCentralityRange  (CentralityMin,CentralityMax);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetDCAparameters    (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetTrackCuts        (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,MaxGoldenChi2,MaxTPCchi2,
                                                                      MaxITSchi2,MaxFracSharedCls,ITSreq);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetPIDCuts          (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetPrefilterCuts    (massLim,phivLim);
    task -> AliAnalysisTaskDielectronsPbPb_Data::SetKinematicCuts    (ptMin,ptMax,etaLim);
    task -> AliAnalysisTaskDielectronsPbPb_Data::EventMixingSettings (MaxNumberEvts,MaxNumberTrks,NumberEvtsToMix,NcentralityBins,NvertexBins,NeventPlaneBins);
    mgr -> AddTask(task);
    AliAnalysisDataContainer *output = mgr->CreateContainer (InputBox,TList::Class(),AliAnalysisManager::kOutputContainer,OutputFile);
    mgr -> ConnectInput  (task,0,input);
    mgr -> ConnectOutput (task,1,output);
    
    
    return task;
}
//____________________________________________________________________________________________________________________________________________________________________________________________________________

