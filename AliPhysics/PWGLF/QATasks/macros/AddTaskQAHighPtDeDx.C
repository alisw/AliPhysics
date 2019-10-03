AliAnalysisTask* AddTaskQAHighPtDeDx(Bool_t AnalysisMC = kFALSE,
			                         Int_t typerun =1, // 0 for pp and 1 for Pb-Pb or pPb
			                         TString  type ="ESD",
			                         UInt_t kTriggerInt = AliVEvent::kINT7, //for pPb kINT7, for pp or PbPb kMB
			                         Float_t minCent = 0.,
                                     Float_t maxCent = 80.,
			                         char *centralityEstimator = "V0A",//for pPb V0A for PbPb V0M
			                         Bool_t ispileuprej = kFALSE
                                    )
{

  /////////////////////////////////////////
  // Few notes to set the wagon in a train
  // For pp collisions:
  //       Bool_t AnalysisMC = kFALSE,
  //       Int_t typerun =0,
  //       TString  type ="AOD",
  //       UInt_t kTriggerInt = AliVEvent::kMB,
  //       Float_t minCent = 0.,//This does not affect pp
  //       Float_t maxCent = 80.,//This does not affect pp
  //       char *centralityEstimator = "V0A",//for pPb V0A for PbPb V0M. This does not affect pp
  //       Bool_t ispileuprej = kFALSE
  // For PbPb collisions:
  //       Bool_t AnalysisMC = kFALSE,
  //       Int_t typerun =1,
  //       TString  type ="AOD",
  //       UInt_t kTriggerInt = AliVEvent::kMB,
  //       Float_t minCent = 0.,
  //       Float_t maxCent = 80.,
  //       char *centralityEstimator = "V0M",
  //       Bool_t ispileuprej = kFALSE
  /////////////////////////////////////////
    
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
  


  
  AliAnalysisTaskQAHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskQAHighPtDeDx("taskHighPtDeDxpp");
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  taskHighPtDeDx->SetAnalysisType(type);
  taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
  if(typerun==1){
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE);
    taskHighPtDeDx->SetMinCent(minCent);
    taskHighPtDeDx->SetMaxCent(maxCent);
    taskHighPtDeDx->SetCentralityEstimator(centralityEstimator);
  }
  else
    taskHighPtDeDx->SetAnalysisPbPb(kFALSE);
  taskHighPtDeDx->SetDebugLevel(0);
  taskHighPtDeDx->SetEtaCut(0.8);
  taskHighPtDeDx->SetVtxCut(10.0);
  taskHighPtDeDx->SetTrigger(kTriggerInt);
  taskHighPtDeDx->SetPileUpRej(ispileuprej);
  //Set Filtesr
  taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
  taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
  taskHighPtDeDx->SetStoreMcIn(AnalysisMC);     // def: kFALSE
  
  mgr->AddTask(taskHighPtDeDx);
  
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //=======================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
  AliAnalysisDataContainer *cout_histdedx;
  cout_histdedx=0;
  cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHighPtDeDx, 1, cout_histdedx);
  
  // Return task pointer at the end
  return taskHighPtDeDx;
   
  
  
}
