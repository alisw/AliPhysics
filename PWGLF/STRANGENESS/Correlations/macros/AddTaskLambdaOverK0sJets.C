AliAnalysisTaskLambdaOverK0sJets *AddTaskLambdaOverK0sJets( TString  name      = "LambdaOverK0sRatio", 
							    TString  data      = "PbPb2010", 
							    Float_t  minCen    = 0.,
							    Float_t  maxCen    = 90.,
							    Float_t  ptMinTrig = 8.,
							    Float_t  ptMaxTrig = 20.,
							    Float_t  etaMaxTrig = 0.75,
							    Float_t  checkIDTrig= kFALSE,
							    Float_t  rapMaxV0  = 0.75,
							    Bool_t   sepInjec  = kTRUE,
							    Bool_t   isMC      = kFALSE,
							    Bool_t   usePID    = kTRUE,
							    Bool_t   doQA      = kFALSE){


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckCascade", "No analysis manager to connect to.");
    return NULL;
  }   
  
  Float_t  nSigmaPID = 3.0;
  
  // Create and configure the task
  AliAnalysisTaskLambdaOverK0sJets *task = new AliAnalysisTaskLambdaOverK0sJets(name.Data());
  task->SetData(data);
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

  AliAnalysisDataContainer *coutput2 =  
    mgr->CreateContainer(name+"_ME", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 name+"_ME.root");

  AliAnalysisDataContainer *coutput3 =  
    mgr->CreateContainer(name+"_QA", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 name+"_QA.root");
  
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  
  return task;
}   
