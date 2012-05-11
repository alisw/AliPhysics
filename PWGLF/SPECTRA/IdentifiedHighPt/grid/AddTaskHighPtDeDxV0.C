AliAnalysisTask* AddTask(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun,  UInt_t kTriggerInt[], Float_t minc[], Float_t maxc[] )
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
  
 //
  // Add track filters, with Golden Cuts
  //
  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = 
    AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);

  
  //old cuts without golden cut
  AliAnalysisFilter* trackFilter0 = new AliAnalysisFilter("trackFilter");
  Bool_t clusterCut=0;
  Bool_t selPrimaries=kTRUE;
  AliESDtrackCuts* esdTrackCutsL0 = new AliESDtrackCuts;
 // TPC  
  if(clusterCut == 0)  esdTrackCutsL0->SetMinNClustersTPC(70);
  else if (clusterCut == 1) {
    esdTrackCutsL0->SetMinNCrossedRowsTPC(70);
    esdTrackCutsL0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }
  else {
    AliWarningClass(Form("Wrong value of the clusterCut parameter (%d), using cut on Nclusters",clusterCut));
    esdTrackCutsL0->SetMinNClustersTPC(70);
  }
  esdTrackCutsL0->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsL0->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsL0->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCutsL0->SetRequireITSRefit(kTRUE);
  esdTrackCutsL0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  if(selPrimaries) {
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCutsL0->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  }
  esdTrackCutsL0->SetMaxDCAToVertexZ(2);
  esdTrackCutsL0->SetDCAToVertex2D(kFALSE);
  esdTrackCutsL0->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsL0->SetMaxChi2PerClusterITS(1e10);
  esdTrackCutsL0->SetMaxChi2TPCConstrainedGlobal(1e10);


  trackFilter0->AddCuts(esdTrackCutsL0);

  
  
  AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
  AliESDtrackCuts* esdTrackCutsTPC = 
    AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackFilterTPC->AddCuts(esdTrackCutsTPC);


  // Create the task and configure it 
  //========================================================================
  if(typerun==2){//heavy ion analysis

    AliAnalysisTaskHighPtDeDxV0* taskHighPtDeDxV0[6];
    for( Int_t i=0; i<6; ++i ){
      taskHighPtDeDxV0[i]=0;
      Char_t TaskName[256]={0};
      sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc[i],maxc[i]);
      
      taskHighPtDeDxV0[i] = new AliAnalysisTaskHighPtDeDxV0(TaskName);
      TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
      
      taskHighPtDeDxV0[i]->SetAnalysisType(type);
      taskHighPtDeDxV0[i]->SetAnalysisMC(AnalysisMC);
      taskHighPtDeDxV0[i]->SetAnalysisPbPb(kTRUE);
      taskHighPtDeDxV0[i]->SetProduceTPCBranch(kFALSE);
      taskHighPtDeDxV0[i]->SetDebugLevel(0);
      taskHighPtDeDxV0[i]->SetEtaCut(0.8);
      taskHighPtDeDxV0[i]->SetVtxCut(10.0);
      taskHighPtDeDxV0[i]->SetMinPt(0.0);           // def: 2.0
      taskHighPtDeDxV0[i]->SetMassCut(0.1);         // def: 0.03
      taskHighPtDeDxV0[i]->SetTreeOption(1);
      taskHighPtDeDxV0[i]->SetRequireRecV0(kFALSE); // def: kTRUE
      taskHighPtDeDxV0[i]->SetStoreMcIn(kTRUE);     // def: kFALSE
      taskHighPtDeDxV0[i]->SetTrigger1(kTriggerInt[0]);
      taskHighPtDeDxV0[i]->SetTrigger2(kTriggerInt[1]);
      taskHighPtDeDxV0[i]->SetMinCent(minc[i]);
      taskHighPtDeDxV0[i]->SetMaxCent(maxc[i]);
      //Set Filtesr
      taskHighPtDeDxV0[i]->SetTrackFilterGolden(trackFilterGolden);
      taskHighPtDeDxV0[i]->SetTrackFilter(trackFilter0);
      taskHighPtDeDxV0[i]->SetTrackFilterTPC(trackFilterTPC);
      
      mgr->AddTask(taskHighPtDeDxV0[i]);
      
      
    }
    
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer *cout_hist[6];
    for( Int_t i=0; i<6; ++i ){
      
      cout_hist[i]=0;
      Char_t outFileName[256]={0};
      sprintf(outFileName,"%s_Tree_%1.0f_%1.0f.root",taskname,minc[i],maxc[i]);
      //AliAnalysisDataContainer *cout_hist    = 0;
      
      cout_hist[i] = mgr->CreateContainer(Form("output_%1.0f_%1.0f",minc[i],maxc[i]), TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
      mgr->ConnectInput (taskHighPtDeDxV0[i], 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskHighPtDeDxV0[i], 1, cout_hist[i]);
      
    }  
    
    // Return task pointer at the end
    return taskHighPtDeDxV0[0];
  }
  
  if(typerun==3){//pp analysis

    AliAnalysisTaskHighPtDeDxV0* taskHighPtDeDxV0 = new AliAnalysisTaskHighPtDeDxV0("taskHighPtDeDxV0pp");
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    taskHighPtDeDxV0->SetAnalysisType(type);
    taskHighPtDeDxV0->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDxV0->SetAnalysisPbPb(kFALSE);
    taskHighPtDeDxV0->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDxV0->SetDebugLevel(0);
    taskHighPtDeDxV0->SetEtaCut(0.8);
    taskHighPtDeDxV0->SetVtxCut(10.0);
    taskHighPtDeDxV0->SetMinPt(0.0);           // def: 2.0
    taskHighPtDeDxV0->SetMassCut(0.03);         // def: 0.03
    taskHighPtDeDxV0->SetTreeOption(1);
    taskHighPtDeDxV0->SetRequireRecV0(kTRUE); // def: kTRUE
    taskHighPtDeDxV0->SetStoreMcIn(kFALSE);     // def: kFALSE
    taskHighPtDeDxV0->SetTrigger1(kTriggerInt[0]);
    taskHighPtDeDxV0->SetTrigger2(kTriggerInt[1]);
    //Set Filtesr
    taskHighPtDeDxV0->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDxV0->SetTrackFilter(trackFilter0);
    taskHighPtDeDxV0->SetTrackFilterTPC(trackFilterTPC);
    
    mgr->AddTask(taskHighPtDeDxV0);
      
   
    
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer *cout_histdedxv0;

    cout_histdedxv0=0;
    Char_t outFileName[256]={0};
    sprintf(outFileName,"%s_Tree.root",taskname);
    //AliAnalysisDataContainer *cout_hist    = 0;
    
    cout_histdedxv0 = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDxV0, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDxV0, 1, cout_histdedxv0);
    
    // Return task pointer at the end
    return taskHighPtDeDxV0;


  }






}
