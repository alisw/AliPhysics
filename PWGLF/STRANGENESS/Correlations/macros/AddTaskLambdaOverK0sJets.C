AliAnalysisTaskLambdaOverK0sJets *AddTaskLambdaOverK0sJets( TString  name      = "LambdaOverK0sRatio", 
							    TString  data      = "PbPb2010", 
							    Float_t  minCen    = 0.,
							    Float_t  maxCen    = 40.,
							    Float_t  fractionSharedTPCcls = 1.,
							    Bool_t   sepInjec  = kTRUE,
							    Bool_t   isMC      = kFALSE,
							    Bool_t   doQA      = kTRUE,
							    Bool_t   useEtaCut = kFALSE){




  Float_t  ptMinTrig   = 5.;
  Float_t  ptMaxTrig   = 1.E100;
  Float_t  etaMaxTrig  = 0.8;
  Float_t  checkIDTrig = kTRUE;
  Float_t  rapMaxV0    = 0.8;
  Bool_t   usePID      = kFALSE;
  Float_t  nSigmaPID   = 3.0;
  Float_t  dcaDaug     = 1.0;
  Float_t  dca2PrmVtx  = 0.1;  
  Float_t  nclsDaug    = 0;
  Float_t  minPtDaughter = 0.;

  Float_t  radiusTPC = 125.;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckCascade", "No analysis manager to connect to.");
    return NULL;
  }   
  
  
  
  // Create and configure the task
  AliAnalysisTaskLambdaOverK0sJets *task = new AliAnalysisTaskLambdaOverK0sJets(name.Data());
  // collision type
  task->SetCollisionType(data);
  task->SetCentrality(minCen,maxCen);
  // trigger particle
  //task->SetTriggerFilterBit(272);
  task->SetTriggerPt(ptMinTrig,ptMaxTrig);
  task->SetTriggerEta(etaMaxTrig);
  task->SetCheckIDTrig(checkIDTrig);
  // V0 candidates
  task->SetEtaCut(useEtaCut);
  task->SetMaxY(rapMaxV0);
  task->SetMaxDCADaughter(dcaDaug); // Added to perform systematics
  task->SetDCAToPrimVtx(dca2PrmVtx); // Added to perform systematics
  //task->SetNSigmaPID(nSigmaPID);
  task->SetNClsTPC(nclsDaug);  // Added to perform systematics
  task->SetMinPtDaughter(minPtDaughter);  
  // PID
  task->SetSeparateInjectedPart(sepInjec);
  // MC
  task->SetMC(isMC);
  task->SetPID(usePID);
  // Setting variables for splitting cut
  task->SetTPCRadius(radiusTPC);    
  task->SetFracSharedTPCcls(fractionSharedTPCcls);
  //task->SetDiffSharedTPCcls(cutSharedTPCcls);
  // QA
  task->SetQA(doQA);
  // Add task
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
                         "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =  
    mgr->CreateContainer(name+"_ME", TList::Class(), 
                         AliAnalysisManager::kOutputContainer, 
                         "AnalysisResults.root");

  AliAnalysisDataContainer *coutput3 =  
    mgr->CreateContainer(name+"_QA", TList::Class(), 
                         AliAnalysisManager::kOutputContainer, 
                         "AnalysisResults.root");
  
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);

  
  return task;
}   
