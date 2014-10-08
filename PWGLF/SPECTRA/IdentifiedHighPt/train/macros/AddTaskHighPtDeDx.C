/*
  Last update: 26/03/2012, vzero branch, in PbPb macro execute the vzero code, a bug was fixed 
  New modifications: 08/10/2014 - more modularity 
*/

AliAnalysisTaskHighPtDeDx* AddTaskHighPtDeDx(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun, Float_t minc, Float_t maxc )
{
  // Creates a pid task and adds it to the analysis manager
  UInt_t kTriggerInt[2] = { AliVEvent::kMB, AliVEvent::kCentral + AliVEvent::kSemiCentral };
    
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
  /*
  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = 
    AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);
  */

  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  esdTrackCutsGolden->SetMinNCrossedRowsTPC(120);
  esdTrackCutsGolden->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsGolden->SetMaxChi2PerClusterITS(36);
  esdTrackCutsGolden->SetMaxFractionSharedTPCClusters(0.4);
  esdTrackCutsGolden->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCutsGolden->SetMaxDCAToVertexXY(3.0);


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
        
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx;
    taskHighPtDeDx=0;
    Char_t TaskName[256]={0};
    sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc,maxc);
     
    taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx(TaskName);
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetMinPt(0.0); // default 2.0
    taskHighPtDeDx->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
    taskHighPtDeDx->SetTreeOption(1);
    taskHighPtDeDx->SetTrigger1(kTriggerInt[0]);
    taskHighPtDeDx->SetTrigger2(kTriggerInt[1]);
    taskHighPtDeDx->SetMinCent(minc);
    taskHighPtDeDx->SetMaxCent(maxc);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
     
    //V0's cuts
    taskHighPtDeDx->SetMassCut(0.1);         // def: 0.03
    taskHighPtDeDx->SetRequireRecV0(kTRUE); // def: kTRUE
    taskHighPtDeDx->SetStoreMcIn(kTRUE);     // def: kFALSE

    mgr->AddTask(taskHighPtDeDx);
    
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer *cout_hist;
    cout_hist=0;
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    //Arrangements for subdirectories (not implemented yet) 
    //outputFileName += ":PWGLF_StrVsMult";
    //if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
     
    cout_hist = mgr->CreateContainer(Form("output_%1.0f_%1.0f",minc,maxc), TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDx, 1, cout_hist);
    
    // Return task pointer at the end
    return taskHighPtDeDx;
  }
}
