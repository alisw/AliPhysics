// $Id$

AliEmcalCompatTask* AddTaskEmcalCompat()
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCompat", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalCompat", "This task requires an input event handler");
    return NULL;
  }
  if (!(mgr->GetInputEventHandler()->InheritsFrom("AliESDInputHandler")))
  {
    ::Error("AddTaskEmcalCompat", "This task works only for (skimmed) ESD");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalCompatTask *ectask = new AliEmcalCompatTask("EmcCompatTask");

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(ectask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  
  mgr->ConnectInput  (ectask, 0,  cinput1 );
  
  return ectask;
}
