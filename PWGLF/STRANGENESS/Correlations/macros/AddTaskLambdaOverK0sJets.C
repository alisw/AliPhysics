AliAnalysisTaskLambdaOverK0sJets *AddTaskLambdaOverK0sJets( TString  name      = "LambdaOverK0sRatio", 
							    Double_t minCen    = 0.,
							    Double_t maxCen    = 90.,
							    Double_t ptMinTrig = 8.,
							    Double_t ptMaxTrig = 20.,
							    Double_t etaMaxTrig = 0.75,
							    Double_t checkIDTrig= kFALSE,
							    Double_t rapMaxV0  = 0.75,
							    Double_t nSigmaPID = 3.0,
							    Bool_t   sepInjec  = kTRUE,
							    Bool_t   isMC      = kFALSE,
							    Bool_t   usePID    = kTRUE,
							    Bool_t   doQA      = kFALSE){


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckCascade", "No analysis manager to connect to.");
    return NULL;
  }   
  
  
  // Create and configure the task
  AliAnalysisTaskLambdaOverK0sJets *task = new AliAnalysisTaskLambdaOverK0sJets(name.Data());
  task->SetCentrality(minCen,maxCen);
  task->SetTriggerPt(ptMinTrig,ptMaxTrig);
  task->SetTriggerEta(etaMaxTrig);
  task->SetCheckIDTrig(checkIDTrig);
  task->SetMaxY(rapMaxV0);
  task->SetNSigmaPID(nSigmaPID);
  task->SetSeparateInjectedPart(sepInjec);
  task->SetMC(isMC);
  task->SetPID(usePID);
  task->SetQA(doQA);
  mgr->AddTask(task);
  
  
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  if (isMC) name+="_mc";
  
  AliAnalysisDataContainer *coutput1 =  
    mgr->CreateContainer(name, TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 name+".root");

  name+="_QA";
  AliAnalysisDataContainer *coutput2 =  
    mgr->CreateContainer(name, TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 name+".root");
  
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  
  return task;
}   
