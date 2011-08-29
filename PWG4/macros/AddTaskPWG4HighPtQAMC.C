void AddTaskPWG4HighPtQAMCAll(char *prodType = "LHC10e14") {

  AliPWG4HighPtQAMC *taskQAMC00 = AddTaskPWG4HighPtQAMC(prodType,0,0);
  AliPWG4HighPtQAMC *taskQAMC00 = AddTaskPWG4HighPtQAMC(prodType,0,1);
  AliPWG4HighPtQAMC *taskQAMC00 = AddTaskPWG4HighPtQAMC(prodType,0,2);
  //  AliPWG4HighPtQAMC *taskQAMC10 = AddTaskPWG4HighPtQAMC(prodType,1,0);
  //  AliPWG4HighPtQAMC *taskQAMC20 = AddTaskPWG4HighPtQAMC(prodType,2,0);
  AliPWG4HighPtQAMC *taskQAMC70 = AddTaskPWG4HighPtQAMC(prodType,7,0);
  AliPWG4HighPtQAMC *taskQAMC70 = AddTaskPWG4HighPtQAMC(prodType,7,1);
}

AliPWG4HighPtQAMC* AddTaskPWG4HighPtQAMC(char *prodType = "LHC10e14", Int_t trackType = 0, Int_t cuts =0)
{
  // Creates HighPtQAMC analysis task and adds it to the analysis manager.

  //Track types:
  //trackType
  //           0: global tracks
  //           1: TPConly tracks
  //           2: TPConly constrained tracks

  //Load common track cut class
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C");
  
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
  AliESDtrackCuts *trackCutsReject = 0x0;
  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    // tight global tracks - RAA analysis
    trackCuts = CreateTrackCutsPWG4(1000);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
    trackCuts = CreateTrackCutsPWG4(10001004);
   }
  if(trackType==0 && cuts==2) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = CreateTrackCutsPWG4(10011004);
  }
  if(trackType==7 && cuts==0) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWG4(10041004);
    trackCutsReject = CreateTrackCutsPWG4(1004);
  }
  if(trackType==7 && cuts==1) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWG4(10011004);
  }
  //Set track cuts for TPConly tracks
  if(trackType==1 || trackType==2) { 
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWG4(2001);
  }
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  
  //Create the task
  AliPWG4HighPtQAMC *taskPWG4QAMC = new AliPWG4HighPtQAMC(Form("AliPWG4HighPtQAMC%d%d",trackType,cuts));
  taskPWG4QAMC->SetCuts(trackCuts);
  taskPWG4QAMC->SetCutsReject(trackCutsReject);
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
  outputfile += Form(":PWG4_HighPtQAMC%d%d",trackType,cuts);
  
  AliAnalysisDataContainer *cout_hist1 = mgr->CreateContainer(Form("qa_histsMC%d%d",trackType,cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);

  mgr->AddTask(taskPWG4QAMC);
  mgr->ConnectInput(taskPWG4QAMC,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4QAMC,0,cout_hist1);

  // Return task pointer at the end
  return taskPWG4QAMC;
}
