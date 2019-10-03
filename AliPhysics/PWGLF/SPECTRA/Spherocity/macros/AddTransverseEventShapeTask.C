AliAnaTransverseEventShapeTask* AddTask(Bool_t AnalysisMC = kFALSE,
					Int_t typerun =0, // 0 for pp and 1 for Pb-Pb or pPb
					TString  type ="ESD",
					UInt_t kTriggerInt = AliVEvent::kMB, //for pPb kINT7, for pp or PbPb kMB
					Float_t minCent = 0.,
					Float_t maxCent = 80.,
					char *centralityEstimator = "V0A",//for pPb V0A for PbPb V0M
					Bool_t ispileuprej = kFALSE
					)
{
 
  
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
  
  gROOT->LoadMacro("$(ALICE_PHYSICS)/PWGJE/macros/CreateTrackCutsPWGJE.C");

  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);

  //Hybrid track cuts, https://twiki.cern.ch/twiki/bin/viewauth/ALICE/HybridTracks
  AliAnalysisFilter* trackFilterHybrid1 = new AliAnalysisFilter("trackFilterHybrid");
  AliESDtrackCuts*  trackCutsHybrid1 = CreateTrackCutsPWGJE(10001006);
  trackFilterHybrid1->AddCuts(trackCutsHybrid1);
  //second part, these tracks need to be constrained
  AliAnalysisFilter* trackFilterHybrid2 = new AliAnalysisFilter("trackFilterHybrid2");
  AliESDtrackCuts*  trackCutsHybrid2 = CreateTrackCutsPWGJE(10041006);
  trackFilterHybrid2->AddCuts(trackCutsHybrid2);
  
  AliAnaTransverseEventShapeTask* taskESA = new AliAnaTransverseEventShapeTask("taskESA");
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  taskESA->SetAnalysisType(type);
  taskESA->SetAnalysisMC(AnalysisMC);
  if(typerun==1){
    taskESA->SetAnalysisPbPb(kTRUE);
    taskESA->SetMinCent(minCent);
    taskESA->SetMaxCent(maxCent);
    taskESA->SetCentralityEstimator(centralityEstimator);
  }
  else
    taskESA->SetAnalysisPbPb(kFALSE);

  //event shape
  taskESA->SetUseHybridESA(kTRUE);
  taskESA->SetTrackFilterESAHyb1(trackFilterHybrid1);
  taskESA->SetTrackFilterESAHyb2(trackFilterHybrid2);
  taskESA->SetTrackFilterESA(trackFilterGolden);
  taskESA->SetMinMultForESA(3);
  taskESA->SetStepSizeESA(0.1);
  taskESA->SetIsEtaAbsESA(kFALSE);
  taskESA->SetTrackEtaMinESA(0.0);
  taskESA->SetTrackEtaMaxESA(0.8);
  taskESA->SetTrackPtMinESA(0.15);
  taskESA->SetTrackPtMaxESA(10.0);
  //

  taskESA->SetDebugLevel(0);
  taskESA->SetEtaCut(0.8);
  taskESA->SetVtxCut(10.0);
  taskESA->SetTrigger(kTriggerInt);
  taskESA->SetPileUpRej(ispileuprej);
  taskESA->SetStoreMcIn(AnalysisMC);     // def: kFALSE

  mgr->AddTask(taskESA);


  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //=======================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cout_histdedx;
  cout_histdedx=0;
  cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput (taskESA, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskESA, 1, cout_histdedx);

  // Return task pointer at the end
  return taskESA;



}
