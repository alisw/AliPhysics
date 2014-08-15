// $Id$

AliEmcalSetupTask* AddTaskEmcalSetup(
  const char *geop    = 0, /*path to geometry folder*/
  const char *oadp    = 0, /*path to OADB folder*/
  const char *ocdp    = 0, /*path to OCDB (if "uselocal", a copy placed in ALIROOT will be used*/
  const char *objs    = 0, /*objects for which alignment should be applied*/
  const Bool_t noOCDB = 0) /*if true then do not mess with OCDB */
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalSetup", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalSetup", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliEmcalSetupTask *eTask = new AliEmcalSetupTask("EmcalSetupTask");
  eTask->SetNoOCDB(noOCDB);
  if (geop) eTask->SetGeoPath(geop);
  if (oadp) eTask->SetOadbPath(oadp);
  if (ocdp) eTask->SetOcdbPath(ocdp);
  if (objs) eTask->SetObjs(objs);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0,  cinput1 );
  return eTask;
}
