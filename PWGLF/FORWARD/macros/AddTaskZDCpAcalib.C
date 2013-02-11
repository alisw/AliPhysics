AliAnalysisTaskSE* AddTaskZDCpAcalib(Bool_t  applyPS = kTRUE,
				  TString outfname = "ZDCpAcalib",
				  Bool_t  isMC = kFALSE)
{
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskZDCPbPb", "No analysis manager to connect to.");
    return NULL;
  }  

  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if(!inputHandler){
    ::Error("AddTaskZDCPbPb", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = inputHandler->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
   AliAnalysisTaskZDCpAcalib* task = new AliAnalysisTaskZDCpAcalib("taskZDCpA");

   if(inputDataType.CompareTo("ESD")==0){
      task->SetInput(1);
      //
      // apply physics selection
      if(applyPS) task->SelectCollisionCandidates();
   }
   else if(inputDataType.CompareTo("AOD")==0){
      task->SetInput(2);
      //printf("  AliAnalysisTaskZDCpAcalib initialized for AOD analysis\n");
   }
   
   if(isMC==kTRUE) task->SetMCInput();

   mgr->AddTask(task);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   AliAnalysisDataContainer *coutput  = mgr->CreateContainer(outfname.Data(), 
   					TList::Class(),
					AliAnalysisManager::kOutputContainer, 	
					Form("%s:ZDCtree", mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, coutput);

   return task;   
}


