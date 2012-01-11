AliAnalysisTaskSE* AddTaskFBFqa(Char_t *name, Bool_t debug=kFALSE)
{
  // Creates a QA task to explore main observables related to FLOW and BALANCE FUNCTION analysis
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFBFqa", "No analysis manager to connect to.");
    return NULL;
  }  
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  // Configure analysis
  //===========================================================================
  //gROOT->LoadMacro("AliGlobalFBFqa.cxx++");
  AliGlobalFBFqa *taskChk = new AliGlobalFBFqa( name );
  if(debug)
    taskChkMB->SetDebugON();
  mgr->ConnectInput (taskChk,0,cinput1);
  AliAnalysisDataContainer *outListMB = mgr->CreateContainer(name,TList::Class(),
                                       AliAnalysisManager::kOutputContainer, 
                                       Form("%s:%s", mgr->GetCommonFileName(),name) );
  mgr->ConnectOutput(taskChk,1,outListMB);

  return taskChk;
}


