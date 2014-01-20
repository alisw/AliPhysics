// $Id$

AliEmcalDebugTask* AddTaskEmcalDebug(
  const UInt_t id       = 0,
  const char *fnametest = 0,
  const char *tnamebase = "DebugTask",
  const char *outputfn  = "AnalysisResults.root"
)
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

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalDebugTask *dtask = new AliEmcalDebugTask(tnamebase);
  UInt_t tid = id;
  if (id==0) {
    TRandom3 r(0);
    tid = r.Integer(kMaxUInt);
  }
  ::Info("AddTaskEmcalDebug","Setting up debug task with id %u", tid);
  dtask->SetId(tid);
  if (fnametest) 
    dtask->SetFileTest(fnametest);
  dtask->SelectCollisionCandidates(0);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(dtask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s_%u",tnamebase,tid),
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputfn);
  mgr->ConnectInput(dtask, 0, cinput);
  mgr->ConnectOutput(dtask,1,coutput);
  
  return dtask;
}
