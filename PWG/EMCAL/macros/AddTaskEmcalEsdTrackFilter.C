// $Id$

AliEmcalEsdTrackFilterTask* AddTaskEmcalEsdTrackFilter(
  const char *name              = "TrackFilter",
  const char *trackCuts         = "Hybrid_LHC11h",
  const char *taskName          = "AliEmcalEsdTrackFilterTask"
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
  strTrackCuts.ToLower();

  if (strTrackCuts.Contains("hybrid")) {
    cutsType = kHybrid;
  } else if (strTrackCuts.Contains("tpconly")) {
    cutsType = kTpcOnly;
    cutsLabel = "TPC only constrained tracks";
  } else {
    ::Warning("AddTaskEmcalEsdTpcTrack", "Cuts type not recognized, will assume Hybrid");
  }

  if (strTrackCuts.Contains("lhc10h")) {
    dataSet = kLHC10h;
  } else if (strTrackCuts.Contains("lhc11a") || strTrackCuts.Contains("lhc12a15a")) {
    dataSet = kLHC11a;
    dataSetLabel = "LHC11a";
  } else if (strTrackCuts.Contains("lhc10e")   ||  strTrackCuts.Contains("lhc10d")   ||
	     strTrackCuts.Contains("lhc10e20") ||  strTrackCuts.Contains("lhc10f6a") ||
	     strTrackCuts.Contains("lhc11a1a") ||  strTrackCuts.Contains("lhc11a1b") ||
	     strTrackCuts.Contains("lhc11a1c") ||  strTrackCuts.Contains("lhc11a1d") ||
	     strTrackCuts.Contains("lhc11a1e") ||  strTrackCuts.Contains("lhc11a1f") ||
	     strTrackCuts.Contains("lhc11a1g") ||  strTrackCuts.Contains("lhc11a1h") ||
	     strTrackCuts.Contains("lhc11a1i") ||  strTrackCuts.Contains("lhc11a1j")) {
    dataSet = kLHC11a;
    dataSetLabel = "LHC10e";
  } else if (strTrackCuts.Contains("lhc11c")) {
    dataSet = kLHC11c;
    dataSetLabel = "LHC11c";
  } else if (strTrackCuts.Contains("lhc11d")) {
    dataSet = kLHC11d;
    dataSetLabel = "LHC11d";
  } else if (strTrackCuts.Contains("lhc11h") || strTrackCuts.Contains("lhc12a15e")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC11h";
  } else if (strTrackCuts.Contains("lhc12g")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC12g";
  } else if (strTrackCuts.Contains("lhc13b")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13b";
  } else if (strTrackCuts.Contains("lhc13c")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13c";
  } else if (strTrackCuts.Contains("lhc13d")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13d";
  } else if (strTrackCuts.Contains("lhc13e")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13e";
  } else if (strTrackCuts.Contains("lhc13f")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13f";
  } else if (strTrackCuts.Contains("lhc13g")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13g";
  } else if (strTrackCuts.Contains("lhc12a15f")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC12a15f";
  } else if (strTrackCuts.Contains("lhc13b4")) {
    dataSet = kLHC11h;
    dataSetLabel = "LHC13b4";
  } else if (strTrackCuts.Contains("lhc12a15g")) {
    dataSet = kLHC11d;
    dataSetLabel = "LHC12a15g";
  } else if (strTrackCuts.Contains("lhc12f2a")) {
    dataSet = kLHC11d;
    dataSetLabel = "LHC12f2a";
  } else {
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
  AliEmcalEsdTrackFilterTask *eTask = new AliEmcalEsdTrackFilterTask(taskName); // default is no cut
  Bool_t includeNoITS = kFALSE;
  if ((dataSet == kLHC11c && cutsType == kHybrid) ||
      (dataSet == kLHC11d && cutsType == kHybrid) ||
      (dataSet == kLHC11h && cutsType == kHybrid)) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001008);       //1000 adds SPD any requirement
    eTask->SetTrackCuts(cutsp);
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041008);       //1004 removes ITSrefit requirement from standard set   
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetHybridTrackCuts(hybsp);
  } else if ((dataSet == kLHC11h && cutsType == kHybrid) ||
	     (dataSet == kLHC11a && cutsType == kHybrid)) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001006);       //1000 adds SPD any requirement
    eTask->SetTrackCuts(cutsp);
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041006);       //1004 removes ITSrefit requirement from standard set    
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetHybridTrackCuts(hybsp);
    includeNoITS = kTRUE;
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
