// $Id$

AliJetTriggerSelectionTask* AddTaskJetTriggerSelection(
  const char *nclusters          = "CaloClusters",
  TF1        *eth                = 0,
  Double_t    maxdistance        = 0.15,
  UInt_t      triggerbits        = AliVEvent::kEMCEJE,
  const char *taskname           = "AliJetTriggerSelectionTask"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetTriggerSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetTriggerSelection", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s",taskname, nclusters));

  AliJetTriggerSelectionTask* jetTask = new AliJetTriggerSelectionTask(name);
  jetTask->AddClusterContainer(nclusters);
  jetTask->SetEnergyThreshold(eth);  
  jetTask->SetMaxDistance(maxdistance);
  jetTask->SetTriggerBits(triggerbits);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  
  return jetTask;
}
