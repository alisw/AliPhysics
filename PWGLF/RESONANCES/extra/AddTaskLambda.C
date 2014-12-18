AliAnalysisTaskLambdaStar *AddTaskLambda
(
      ULong64_t Triggermask = AliVEvent::kCentral,
      Bool_t Cirpid = kFALSE,
      Int_t Centmin = 0,
      Int_t Centmax =10,
      Double_t Nsigma = 3.0,
      Int_t Nmix = 2,
      Bool_t Nstrong = kFALSE,
      Int_t Centp = 510
)
 {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;}  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL; }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "AOD") {
    Printf("ERROR! This task can only run on AODs!");}

  // Configure analysis
   AliAnalysisTaskLambdaStar *task = new AliAnalysisTaskLambdaStar("LStar");
  task->SelectCollisionCandidates(Triggermask);
  task->SetCentrality(Centmin,Centmax );
  task->UseCircPID(Cirpid);
  task->SetNMix(Nmix);
  task->SetNSigma(Nsigma);
  task->SetStrongCut(Nstrong);
  task->SetCentPatch(Centp);
   
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cHistLambda", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:lambdastar", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cHistPrim", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:primary", AliAnalysisManager::GetCommonFileName()));
    
	mgr->ConnectInput (task, 0, cinput0);
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);
    
  return task;
}   


