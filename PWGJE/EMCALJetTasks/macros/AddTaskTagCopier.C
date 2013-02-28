// $Id: AddTaskSAQA.C 60163 2013-01-03 09:37:04Z loizides $

AliJetConstituentTagCopier* AddTaskTagCopier(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *nmcparticles       = "MCParticles",
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
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
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "_TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "_EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "_USER";

  AliJetConstituentTagCopier* task = new AliJetConstituentTagCopier(name);
  task->SetTracksName(ntracks);
  task->SetClusName(nclusters);
  task->SetMCParticlesName(nmcparticles);
  task->SetTrackPtCut(trackptcut);
  task->SetClusPtCut(clusptcut);
  task->SetAnaType(type);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput(task, 0, cinput1);

  return task;
} 
