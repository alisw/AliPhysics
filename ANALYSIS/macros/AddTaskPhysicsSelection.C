AliPhysicsSelectionTask* AddTaskPhysicsSelection(Bool_t mCAnalysisFlag = kFALSE, Bool_t deprecatedFlag = kTRUE, UInt_t computeBG = 0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Configure analysis
  //===========================================================================
    
    

  AliPhysicsSelectionTask *task = new AliPhysicsSelectionTask("");
  mgr->AddTask(task);
  
  AliPhysicsSelection* physSel = task->GetPhysicsSelection();
  if (mCAnalysisFlag)      
    physSel->SetAnalyzeMC();
  if (computeBG)
    physSel->SetComputeBG(computeBG);

  if(!deprecatedFlag) 
    AliFatal("The BG ID flag is deprecated. Please use the OADB to configure the cuts");

  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cstatsout",
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
		
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);

  mgr->RegisterExtraFile("event_stat.root");

  return task;
}   
