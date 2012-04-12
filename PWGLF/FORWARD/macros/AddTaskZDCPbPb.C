AliAnalysisTaskSE* AddTaskZDCPbPb(Float_t centrlowlim = 0.,
                                  Float_t centruplim = 90.,
                                  TString centrest = "V0M",
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
   gROOT->LoadMacro("AliAnalysisTaskZDCPbPb.cxx++g");   
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
   
   if(isMC==kTRUE) task->SetMCInput();
   
   // apply physics selection
   //task->SelectCollisionCandidates();

   mgr->AddTask(task);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //outputFileName += ":outputZDCPbPb";
   
   Printf("AddTaskZDCPbPb - Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutput  = mgr->CreateContainer("ZDChistos", TList::Class(),
					AliAnalysisManager::kOutputContainer, 								Form("%s:ZDCPbPb", mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, coutput);

   return task;   
}


