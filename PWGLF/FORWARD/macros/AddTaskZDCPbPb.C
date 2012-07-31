AliAnalysisTaskSE* AddTaskZDCPbPb(Bool_t  applyPS = kFALSE,
				  Float_t centrlowlim = 0.,
                                  Float_t centruplim = 100.,
                                  TString centrest = "V0M",
				  TString outfname = "ZDCPbPb",
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
   AliAnalysisTaskZDCPbPb* task = new AliAnalysisTaskZDCPbPb("taskZDCPbPb");

   if(inputDataType.CompareTo("ESD")==0){
      task->SetInput(1);
      printf("  AliAnalysisTaskZDCPbPb initialized for ESD analysis\n");
   }
   else if(inputDataType.CompareTo("AOD")==0){
      task->SetInput(2);
      printf("  AliAnalysisTaskZDCPbPb initialized for AOD analysis\n");
   }
   task->SetCentralityRange(centrlowlim, centruplim);
   task->SetCentralityEstimator(centrest);

   
   // apply physics selection
   if(applyPS) task->SelectCollisionCandidates();
   
   if(isMC==kTRUE) task->SetMCInput();

   mgr->AddTask(task);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   AliAnalysisDataContainer *coutput  = mgr->CreateContainer(outfname.Data(), 
   					TList::Class(),
					AliAnalysisManager::kOutputContainer, 	
					Form("%s:ZDCHistos", mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, coutput);

   return task;   
}


