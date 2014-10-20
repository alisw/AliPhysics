// $Id$

AliAnalysisTaskEmcalSample* AddTaskEmcalSample(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  Int_t       nCentBins          = 1,
  const char *taskname           = "AliAnalysisTaskEmcalSample"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalSample", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalSample", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalSample* emcTask = new AliAnalysisTaskEmcalSample(name);
  emcTask->SetCentRange(0.,100.);
  emcTask->SetNCentBins(nCentBins);

  AliParticleContainer *trackCont  = emcTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = emcTask->AddClusterContainer(nclusters);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(emcTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (emcTask, 0,  cinput1 );
  mgr->ConnectOutput (emcTask, 1, coutput1 );
  
  return emcTask;
}
