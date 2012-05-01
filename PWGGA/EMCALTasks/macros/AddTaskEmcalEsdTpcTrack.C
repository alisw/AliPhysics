// $Id$

AliEmcalEsdTpcTrackTask* AddTaskEmcalEsdTpcTrack(
  const char *name       = "TpcSpdVertexConstrainedTracks"
)
{  
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

  AliESDtrackCuts *cutsp = new AliESDtrackCuts;
  /*hybrid tracks*/
  // TPC
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  cutsp->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  cutsp->SetMinNClustersTPC(70);
  cutsp->SetMaxChi2PerClusterTPC(4);
  cutsp->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  cutsp->SetAcceptKinkDaughters(kFALSE);
  cutsp->SetRequireTPCRefit(kTRUE);
  cutsp->SetMaxFractionSharedTPCClusters(0.4);
  // ITS
  cutsp->SetRequireITSRefit(kTRUE);
  //accept secondaries
  cutsp->SetMaxDCAToVertexXY(2.4);
  cutsp->SetMaxDCAToVertexZ(3.2);
  cutsp->SetDCAToVertex2D(kTRUE);
  //reject fakes
  cutsp->SetMaxChi2PerClusterITS(36);
  cutsp->SetMaxChi2TPCConstrainedGlobal(36);
  cutsp->SetRequireSigmaToVertex(kFALSE);
  cutsp->SetEtaRange(-0.9,0.9);
  cutsp->SetPtRange(0.15, 1E+15);
  //tag = "Global tracks jet analysis with ITSrefit and NclsIter1=PtDep, noSPD requirement, no upper pt cut, golden chi2";
  hybsp = new AliESDtrackCuts(*cutsp);
  cutsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  //tag += " + additonal: SPD any requirement";
  hybsp->SetRequireITSRefit(kFALSE);
  //tag += " + additional: ITSrefit=kFALSE";

  AliEmcalEsdTpcTrackTask *eTask = new AliEmcalEsdTpcTrackTask();
  eTask->SetTrackCuts(cutsp);
  eTask->SetHybridTrackCuts(hybsp);
  eTask->SetTracksName(name);

  cout << " *** TPC track to SPD vertex task configured *** " << endl;

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
}
