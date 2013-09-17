// $Id$

AliJetConstituentTagCopier* AddTaskTagCopier(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *nmcparticles       = "MCParticles",
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  const char *taskname           = "AliJetConstituentTagCopier"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskTagCopier", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskTagCopier", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  
  TString name(taskname);
  if (strcmp(ntracks,"")) {
    name += "_";
    name += ntracks;
  }
  if (strcmp(nclusters,"")) {
    name += "_";
    name += nclusters;
  }
  if (strcmp(nmcparticles,"")) {
    name += "_";
    name += nmcparticles;
  }

  AliJetConstituentTagCopier* task = new AliJetConstituentTagCopier(name);

  AliParticleContainer *trackCont = jetTask->AddParticleContainer(ntracks);
  if (trackCont) trackCont->SetParticlePtCut(trackptcut);

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);
  if (clusterCont) clusterCont->SetClusPtCut(clusptcut);

  AliParticleContainer *mcPartCont = jetTask->AddParticleContainer(nmcparticles);

  task->ConnectMCParticleContainerID(mcPartCont);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput(task, 0, cinput1);

  return task;
} 
