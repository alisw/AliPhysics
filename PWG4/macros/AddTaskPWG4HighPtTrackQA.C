void AddTaskPWG4HighPtTrackQAAll(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  int cent = 10;
  
  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA02cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
  // AliPWG4HighPtTrackQA *taskTrackQA10cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
  // AliPWG4HighPtTrackQA *taskTrackQA11cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA20cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA21cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA40cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA41cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA50cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,5,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA60cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,6,0);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
      //    AliPWG4HighPtTrackQA *taskTrackQA02 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
      // AliPWG4HighPtTrackQA *taskTrackQA10 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
      // AliPWG4HighPtTrackQA *taskTrackQA11 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA20 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA21 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA40 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA41 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA50 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,5,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA60 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,6,0);
      AliPWG4HighPtTrackQA *taskTrackQA70 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);
    }
  }

}

void AddTaskPWG4HighPtTrackQAAllReduced(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  int cent = 10;
  
  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
      AliPWG4HighPtTrackQA *taskTrackQA70 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);
    }
  }

}

void AddTaskPWG4HighPtTrackQAAOD(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 1, Int_t filterBit) 
{   
  AliPWG4HighPtTrackQA *taskTrackQA = AddTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,0,0,0);
  taskTrackQA->SetFilterMask(filterBit);
}

AliPWG4HighPtTrackQA* AddTaskPWG4HighPtTrackQA(char *prodType = "LHC10e14",Bool_t isPbPb=kTRUE,Int_t iAODanalysis = 0, Int_t centClass = 0, Int_t trackType = 0, Int_t cuts = 0)
{
  /*
    trackType: 0 = global
               1 = TPC stand alone
               2 = TPC stand alone constrained to SPD vertex
               4 = TPC stand alone constrained to SPD vertex with QA track selection on global tracks
	       5 = Hybrid tracks: constrained TPConly for which no tight ITS is available
               6 = Hybrid tracks: constrained loose global for which no tight ITS is available
    cuts:      0 (global) = standard ITSTPC2010 a la RAA analysis
               1 (global) = ITSrefit, no SPD requirements -> standard for jet analysis
               2 (global) = ITSrefit + no hits in SPD
	       3 (global) = standard ITS tight cuts with nCrossed rows cut for hybrid tracks
               0 (TPC)    = standard TPC + NClusters>70
               1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations
               0 (hybrid 5) = constrained TPConly for which no tight ITS is available
               0 (hybrid 6) = constrained loose global for which no tight ITS is available
   */

  //Load common track cut class
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C");

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
  AliESDtrackCuts *trackCutsReject = 0x0;
  AliESDtrackCuts *trackCutsTPConly = new AliESDtrackCuts("AliESDtrackCutsTPConly","TPC only Cuts");

  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    // tight global tracks - RAA analysis
    trackCuts = CreateTrackCutsPWG4(1000);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
    trackCuts = CreateTrackCutsPWG4(10001005);
  }
  if(trackType==0 && cuts==2) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = CreateTrackCutsPWG4(10011005);
  }
  if(trackType==7 && cuts==0) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWG4(10041005);
    trackCutsReject = CreateTrackCutsPWG4(1005);
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  if(trackType==7 && cuts==1) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWG4(10011005);
  }
  if(trackType==7 && cuts==2) {
    // no requirements on SPD and ITSrefit failed
    trackCuts = CreateTrackCutsPWG4(10041005);       //no ITSrefit requirement filter 256
    trackCutsReject = CreateTrackCutsPWG4(10001005); //ITSrefit requirement filter 16
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }

  if(trackType==1 && cuts==0) {
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  if(trackType==1 && cuts==1) {
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWG4(10032001);
  }

  if(trackType==2 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  if(trackType==2 && cuts==1) {
    //Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWG4(10032001);
  }

  if(trackType==4 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  if(trackType==4 && cuts==1) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWG4(10032001);
  }
  if(trackType==5 || trackType==6) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWG4(1003);

    trackCutsReject = CreateTrackCutsPWG4(10021003); 
    
    trackCutsTPConly = CreateTrackCutsPWG4(2002);

    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
    
    trackCutsTPConly->SetEtaRange(-0.9,0.9);
    trackCutsTPConly->SetPtRange(0.15, 1e10);

  }

  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  


  //Create the task
  AliPWG4HighPtTrackQA *taskPWG4TrackQA = new AliPWG4HighPtTrackQA(Form("AliPWG4HighPtTrackQACent%dTrack%dCuts%d",centClass,trackType,cuts));
  taskPWG4TrackQA->SetTrackType(trackType);
  taskPWG4TrackQA->SetCuts(trackCuts);
  taskPWG4TrackQA->SetCutsITSLoose(trackCutsReject);
  taskPWG4TrackQA->SetCutsTPConly(trackCutsTPConly);
  
  taskPWG4TrackQA->SetPtMax(100.);
 
  if(iAODanalysis)
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kAOD);
  else
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kESD);

  if(isPbPb) {
    taskPWG4TrackQA->SetIsPbPb(kTRUE);
    taskPWG4TrackQA->SetCentralityClass(centClass);
  }
  //  taskPWG4TrackQA->SetSigmaConstrainedMax(5.);

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
