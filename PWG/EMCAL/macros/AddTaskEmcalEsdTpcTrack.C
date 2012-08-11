// $Id$

AliEmcalEsdTpcTrackTask* AddTaskEmcalEsdTpcTrack(
						 const char *name       = "TpcSpdVertexConstrainedTracks",
						 const char *trackCuts  = "Hybrid_LHC11h"
						 )
{ 
  enum CutsType {
    kHybrid  = 0,
    kTpcOnly = 1
  };

  enum DataSet {
    kLHC11h  = 0,
    kLHC11a  = 1
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

  if (strTrackCuts.Contains("LHC11h")) {
    dataSet = kLHC11h;
  }
  else if (strTrackCuts.Contains("LHC11a")) {
    dataSet = kLHC11a;
    dataSetLabel = "LHC11a";
  }
  else {
    ::Warning("AddTaskEmcalEsdTpcTrack", "Dataset not recognized, will assume LHC11h");
  }

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalEsdTpcTrack", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalEsdTpcTrack", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");

  AliEmcalEsdTpcTrackTask *eTask = new AliEmcalEsdTpcTrackTask(); // default is TPC only tracks constrained to the vertex

  if (dataSet == kLHC11h && cutsType == kHybrid) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001007);       //1000 adds SPD any requirement
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041007);       //1004 removes ITSrefit requirement from standard set    
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetTrackCuts(cutsp);
    eTask->SetHybridTrackCuts(hybsp);
  }
  else if (dataSet == kLHC11a && cutsType == kHybrid) {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001006);       //1000 adds SPD any requirement
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041006);       //1004 removes ITSrefit requirement from standard set    
    hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    eTask->SetTrackCuts(cutsp);
    eTask->SetHybridTrackCuts(hybsp);
  }

  eTask->SetTracksName(name);

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
