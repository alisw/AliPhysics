// $Id$

AliEmcalEsdTpcTrackTask* AddTaskEmcalEsdTpcTrack(
  const char *name              = "TpcSpdVertexConstrainedTracks",
  const char *trackCuts         = "Hybrid_LHC11h",
  Bool_t      includeNoITS      = kFALSE
)
{ 
  enum CutsType {
    kHybrid  = 0,
    kTpcOnly = 1
  };

  enum DataSet {
    kLHC10h  = 0,
    kLHC11a  = 1,
    kLHC11c  = 3,
    kLHC11d  = 3,
    kLHC11h  = 3
  };

  CutsType cutsType = kHybrid;
  DataSet  dataSet  = kLHC11h;

  TString cutsLabel("hybrid tracks");
  TString dataSetLabel("LHC11h");

  TString strTrackCuts(trackCuts);

  if (strTrackCuts.Contains("Hybrid")) {
    cutsType = kHybrid;
  }
  else if (strTrackCuts.Contains("TpcOnly")) {
    cutsType = kTpcOnly;
    cutsLabel = "TPC only constrained tracks";
  }
  else {
    ::Warning("AddTaskEmcalEsdTpcTrack", "Cuts type not recognized, will assume Hybrid");
  }

  if (strTrackCuts.Contains("LHC10h")) {
    dataSet = kLHC10h;
  }
  else if (strTrackCuts.Contains("LHC11a") || strTrackCuts.Contains("LHC12a15a")) {
    dataSet = kLHC11a;
    dataSetLabel = "LHC11a";
  }
  else if (strTrackCuts.Contains("LHC11c")) {
    dataSet = kLHC11c;
    dataSetLabel = "LHC11c";
  }
  else if (strTrackCuts.Contains("LHC11d")) {
    dataSet = kLHC11d;
    dataSetLabel = "LHC11d";
  }
  else if (strTrackCuts.Contains("LHC11h") || strTrackCuts.Contains("LHC12a15e"))
  {
    dataSet = kLHC11h;
    dataSetLabel = "LHC11h";
  }
  else if (strTrackCuts.Contains("LHC12g"))
  {
    dataSet = kLHC11h;
    dataSetLabel = "LHC12g";
  }
  else {
    ::Warning("AddTaskEmcalEsdTpcTrack", "Dataset not recognized, will assume LHC11h");
  }

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalEsdTpcTrack", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalEsdTpcTrack", "This task requires an input event handler");
    return NULL;
  }
  
  if (!evhand->InheritsFrom("AliESDInputHandler")) {
    ::Info("AddTaskEmcalEsdTpcTrack", "This task is only needed for ESD analysis. No task added.");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");

  AliEmcalEsdTpcTrackTask *eTask = new AliEmcalEsdTpcTrackTask(); // default is TPC only tracks constrained to the vertex

  if ((dataSet == kLHC11c && cutsType == kHybrid) ||
      (dataSet == kLHC11d && cutsType == kHybrid) ||
      (dataSet == kLHC11h && cutsType == kHybrid)) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001008);       //1000 adds SPD any requirement
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041008);       //1004 removes ITSrefit requirement from standard set   
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetTrackCuts(cutsp);
    eTask->SetHybridTrackCuts(hybsp);
  }
  else if ((dataSet == kLHC10h && cutsType == kHybrid) ||
           (dataSet == kLHC11a && cutsType == kHybrid)) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001006);       //1000 adds SPD any requirement
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041006);       //1004 removes ITSrefit requirement from standard set    
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetTrackCuts(cutsp);
    eTask->SetHybridTrackCuts(hybsp);
  }

  eTask->SetTracksName(name);
  eTask->SetIncludeNoITS(includeNoITS);

  cout << " *** Track selector task configured to select " << cutsLabel  << " in dataset "<< dataSetLabel << " *** " << endl;

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1);
  
  return eTask;
}
