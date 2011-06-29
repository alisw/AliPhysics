
//DEFINITION OF A FEW CONSTANTS

AliPWG4HighPtQAMC* AddTaskPWG4HighPtQAMC(char *prodType = "LHC10e14", int trackType = 0)
{
  // Creates HighPtQAMC analysis task and adds it to the analysis manager.
  
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
    ::Error("AddPWG4TaskHighPtQAMC", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0) {
    trackCuts = trackCuts->GetStandardITSTPCTrackCuts2010(kTRUE);//Primary Track Selection
    trackCuts->SetRequireITSRefit(kTRUE);
  }
  if(trackType==3) {
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
    
    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1e10);
  }
  //Set track cuts for TPConly tracks
  if(trackType==1 || trackType==2) { 
    trackCuts = trackCuts->GetStandardTPCOnlyTrackCuts(); 
    trackCuts->SetMinNClustersTPC(70);
  }
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  
  AliESDtrackCuts *trackCutsITS = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts with SPD or SDD");
  //Cuts SPD || SDD
  // TPC  
  trackCutsITS->SetMinNClustersTPC(70);
  trackCutsITS->SetMaxChi2PerClusterTPC(4);
  trackCutsITS->SetAcceptKinkDaughters(kFALSE);
  trackCutsITS->SetRequireTPCRefit(kTRUE);
  // ITS
  trackCutsITS->SetRequireITSRefit(kTRUE);
  trackCutsITS->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  trackCutsITS->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);

  trackCutsITS->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  trackCutsITS->SetMaxDCAToVertexZ(2);
  trackCutsITS->SetDCAToVertex2D(kFALSE);
  trackCutsITS->SetRequireSigmaToVertex(kFALSE);

  trackCutsITS->SetEtaRange(-0.9,0.9);
  trackCutsITS->SetPtRange(0.15, 1e10);
  trackCutsITS->SetRequireITSRefit(kTRUE);

  //Create the task
  AliPWG4HighPtQAMC *taskPWG4QAMC = new AliPWG4HighPtQAMC(Form("AliPWG4HighPtQAMC%d",trackType));
  taskPWG4QAMC->SetCuts(trackCuts);
  taskPWG4QAMC->SetCutsITS(trackCutsITS);
  taskPWG4QAMC->SetTrackType(trackType);
  
  if(!strcmp(prodType, "LHC10e14")) taskPWG4QAMC->SetPtMax(500.);
  else taskPWG4QAMC->SetPtMax(100.);

  //taskPWG4QAMC->SetSigmaConstrainedMax(5.);

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  //  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  printf("Create output containers \n");
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtQAMC%d",trackType);
  
  AliAnalysisDataContainer *cout_hist1 = mgr->CreateContainer(Form("qa_histsMC%d",trackType), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *cout_hist2 = mgr->CreateContainer(Form("qa_histsMCITS%d",trackType), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);  

  mgr->AddTask(taskPWG4QAMC);
  mgr->ConnectInput(taskPWG4QAMC,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4QAMC,0,cout_hist1);
  mgr->ConnectOutput(taskPWG4QAMC,1,cout_hist2);

    // Return task pointer at the end
  return taskPWG4QAMC;
}
