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

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  TString inputDataType = inputHandler->GetDataType(); // can be "ESD" or "AOD"

  // Configure analysis
  //===========================================================================
  AliPhysicsSelectionTask *task = new AliPhysicsSelectionTask("");

  // this makes physics selection to work using AliMultiInputEventHandler
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    AliInputEventHandler *ih = multiInputHandler->GetFirstInputEventHandler();
    if (!ih) {
      ::Error("AddTaskPhysicsSelection","ESD or AOD input handler is missing");
      return NULL;
    }
    ih->SetEventSelection(multiInputHandler->GetEventSelection());
    inputDataType = ih->GetDataType(); // can be "ESD" or "AOD"
  }
  
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

  return task;
}   
