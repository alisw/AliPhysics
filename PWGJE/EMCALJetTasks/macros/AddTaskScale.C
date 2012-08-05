// $Id$

AliAnalysisTaskScale* AddTaskScale(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClustersCorr",
  Double_t    ptcut          = 0.150,
  const char *outfilename    = "AnalysisResults.root",
  const char *taskname       = "Scale"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskScale", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskScale", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s_%d", taskname, nTracks, nClusters, TMath::FloorNint(ptcut*1000)));
  AliAnalysisTaskScale *scaletask = new AliAnalysisTaskScale(name);
  scaletask->SetTracksName(nTracks);
  scaletask->SetClusName(nClusters);
  scaletask->SetPtCut(ptcut);
  scaletask->SetAnaType(AliAnalysisTaskEmcal::kEMCAL);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(scaletask);

  // Create containers for input/output
  TString contname(name);
  contname += "_Histos";
  mgr->ConnectInput (scaletask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coscale = mgr->CreateContainer(contname,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(scaletask,1,coscale);

  return scaletask;
}
