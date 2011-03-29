void AddTaskPWG4HighPtTrackQAAll(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    
  int cent = 10;
  
  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
  AliPWG4HighPtTrackQA *taskTrackQA02cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
  AliPWG4HighPtTrackQA *taskTrackQA10cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
  AliPWG4HighPtTrackQA *taskTrackQA11cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
  AliPWG4HighPtTrackQA *taskTrackQA20cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
  AliPWG4HighPtTrackQA *taskTrackQA21cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
    
  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
      AliPWG4HighPtTrackQA *taskTrackQA02 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
      AliPWG4HighPtTrackQA *taskTrackQA10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
      AliPWG4HighPtTrackQA *taskTrackQA11 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
      AliPWG4HighPtTrackQA *taskTrackQA20 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
      AliPWG4HighPtTrackQA *taskTrackQA21 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
    }
  }

}

AliPWG4HighPtTrackQA* AddTaskPWG4HighPtTrackQA(char *prodType = "LHC10e14",Bool_t isPbPb=kTRUE,Int_t iAODanalysis = 0, Int_t centClass = 0, Int_t trackType = 0, Int_t cuts = 0)
{
  /*
    trackType: 0 = global
               1 = TPC stand alone
               2 = TPC stand alone constrained to SPD vertex
    cuts:      0 (global) = standard ITSTPC2010
               1 (global) = ITSrefit, no SPD requirements
               2 (global) = SPD || SDD
               0 (TPC)    = standard TPC + NClusters>70
               1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations
   */

  // Creates HighPtTrackQA analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtQMC", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddPWG4TaskHighPtTrackQA", "This task requires an input event handler");
    return NULL;
  }  

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    trackCuts = trackCuts->GetStandardITSTPCTrackCuts2010(kTRUE);//Primary Track Selection
    trackCuts->SetRequireITSRefit(kTRUE);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement
    // TPC  
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    
    trackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackCuts->SetMaxDCAToVertexZ(2);
    trackCuts->SetDCAToVertex2D(kFALSE);
    trackCuts->SetRequireSigmaToVertex(kFALSE);
  }
  if(trackType==0 && cuts==2) {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts with SPD or SDD");
    //Cuts SPD || SDD
    // TPC  
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
    
    trackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackCuts->SetMaxDCAToVertexZ(2);
    trackCuts->SetDCAToVertex2D(kFALSE);
    trackCuts->SetRequireSigmaToVertex(kFALSE);
    
    trackCuts->SetRequireITSRefit(kTRUE);
  }
  if(trackType==1 && cuts==0) {
    //Set track cuts for TPConly tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(70);
  }
  if(trackType==1 && cuts==1) {
    //Set track cuts for TPConly tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(0);
  }

  if(trackType==2 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts();
    trackCuts->SetMinNClustersTPC(70);
  }
  if(trackType==2 && cuts==1) {
    //Set track cuts for TPConly constrained tracks
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts();
    trackCuts->SetMinNClustersTPC(0);
  }
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  


  //Create the task
  AliPWG4HighPtTrackQA *taskPWG4TrackQA = new AliPWG4HighPtTrackQA(Form("AliPWG4HighPtTrackQACent%dTrack%dCuts%d",centClass,trackType,cuts));
  taskPWG4TrackQA->SetCuts(trackCuts);
  taskPWG4TrackQA->SetTrackType(trackType);
  
  taskPWG4TrackQA->SetPtMax(100.);
 
  if(iAODanalysis)
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kAOD);
  else
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kESD);

  if(isPbPb) {
    taskPWG4TrackQA->SetIsPbPb(kTRUE);
    taskPWG4TrackQA->SetCentralityClass(centClass);
  }

  taskPWG4TrackQA->SelectCollisionCandidates();


  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtTrackQACent%dTrackType%dCuts%d",centClass,trackType,cuts);
  
  AliAnalysisDataContainer *cout_histQAtrack = mgr->CreateContainer(Form("qa_histsQAtrackCent%dType%dcuts%d",centClass,trackType,cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);

  mgr->AddTask(taskPWG4TrackQA);
  mgr->ConnectInput(taskPWG4TrackQA,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4TrackQA,1,cout_histQAtrack);

  // Return task pointer at the end
  return taskPWG4TrackQA;
}
