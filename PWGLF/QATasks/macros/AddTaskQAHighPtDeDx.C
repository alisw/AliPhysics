

AliAnalysisTask* AddTask(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun, Float_t minc[], Float_t maxc[] )
{
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHighPtDeDx", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskHighPtDeDx", "This task requires an input event handler");
    return NULL;
  }  
  



  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);
  
  AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
  AliESDtrackCuts* esdTrackCutsTPC = 
    AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackFilterTPC->AddCuts(esdTrackCutsTPC);
  


  // Create the task and configure it 
  //========================================================================
  if(typerun==2){//heavy ion analysis
    
    
    AliAnalysisTaskQAHighPtDeDx* taskHighPtDeDx[6];
    for( Int_t i=0; i<6; ++i ){
      taskHighPtDeDx[i]=0;
      Char_t TaskName[256]={0};
      sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc[i],maxc[i]);
      
      taskHighPtDeDx[i] = new AliAnalysisTaskQAHighPtDeDx(TaskName);
      TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
      taskHighPtDeDx[i]->SetAnalysisType(type);
      taskHighPtDeDx[i]->SetAnalysisMC(AnalysisMC);
      taskHighPtDeDx[i]->SetAnalysisPbPb(kTRUE);
      taskHighPtDeDx[i]->SetDebugLevel(0);
      taskHighPtDeDx[i]->SetEtaCut(0.8);
      taskHighPtDeDx[i]->SetVtxCut(10.0);
      taskHighPtDeDx[i]->SetCentralityEstimator(centralityEstimator);
      taskHighPtDeDx[i]->SetTrigger(kTriggerInt);
      taskHighPtDeDx[i]->SetMinCent(minc[i]);
      taskHighPtDeDx[i]->SetMaxCent(maxc[i]);
      taskHighPtDeDx[i]->SetPileUpRej(ispileuprej);
  
      taskHighPtDeDx[i]->SetTrackFilterGolden(trackFilterGolden);
      taskHighPtDeDx[i]->SetTrackFilterTPC(trackFilterTPC);
      taskHighPtDeDx[i]->SetStoreMcIn(analysisMC);

      mgr->AddTask(taskHighPtDeDx[i]);
      
    }
    
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer *cout_hist[6];
    for( Int_t i=0; i<6; ++i ){
      
      cout_hist[i]=0;
      Char_t outFileName[256]={0};
      sprintf(outFileName,"%s_%1.0f_%1.0f.root",taskname,minc[i],maxc[i]);
      //AliAnalysisDataContainer *cout_hist    = 0;
      
      cout_hist[i] = mgr->CreateContainer(Form("output_%1.0f_%1.0f",minc[i],maxc[i]), TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
      mgr->ConnectInput (taskHighPtDeDx[i], 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskHighPtDeDx[i], 1, cout_hist[i]);
      
    }  
    
    // Return task pointer at the end
    return taskHighPtDeDx[0];
  }
  if(typerun==3){//pp analysis
  
    
    AliAnalysisTaskQAHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskQAHighPtDeDx("taskHighPtDeDxpp");
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetTrigger(kTriggerInt);
    taskHighPtDeDx->SetPileUpRej(ispileuprej);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    taskHighPtDeDx->SetStoreMcIn(analysisMC);     // def: kFALSE

    mgr->AddTask(taskHighPtDeDx);
      
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
  
    AliAnalysisDataContainer *cout_histdedx;

    cout_histdedx=0;
    Char_t outFileName[256]={0};
    sprintf(outFileName,"%s.root",taskname);
    //AliAnalysisDataContainer *cout_hist    = 0;
 
    cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDx, 1, cout_histdedx);
 
    // Return task pointer at the end
    return taskHighPtDeDx;
    

    
    
    
  }



 if(typerun==4){//pPb analysis
   cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<                 Runing pPb    <<<<<<<<<<<<<<"<<endl;
    
    AliAnalysisTaskQAHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskQAHighPtDeDx("taskHighPtDeDxpp");
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE);
    taskHighPtDeDx->SetMinCent(-2000);
    taskHighPtDeDx->SetMaxCent(2000);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetTrigger(kTriggerInt);
    taskHighPtDeDx->SetPileUpRej(ispileuprej);
    taskHighPtDeDx->SetCentralityEstimator(centralityEstimator);
    taskHighPtDeDx->SetStoreMcIn(analysisMC);     // def: kFALSE
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
 
    mgr->AddTask(taskHighPtDeDx);
      
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
  
    AliAnalysisDataContainer *cout_histdedx;

    cout_histdedx=0;
    Char_t outFileName[256]={0};
    sprintf(outFileName,"%s.root",taskname);
    //AliAnalysisDataContainer *cout_hist    = 0;
 
    cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDx, 1, cout_histdedx);
 
    // Return task pointer at the end
    return taskHighPtDeDx;
    

    
    
    
  }



  
  
}
