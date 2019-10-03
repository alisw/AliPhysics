AliAnalysisTaskLambdaStar *AddTaskLambda
(
      UInt_t Triggermask = AliVEvent::kCentral,
      Bool_t Cirpid = kFALSE,
      Int_t Centmin = 0,
      Int_t Centmax =10,
      Double_t Nsigma = 3.0,
      Double_t Rejnsigma =0.0,
      Int_t Nmix = 5,
      Int_t Centp = 510,
      Int_t ClusterTPC = 70,
      Float_t DCAxy = 0.1,
      Int_t FilterBit= 032,
      TString  Dirsuffixname = ""
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
  task->SetRejNSigma(Rejnsigma);
  task->SetCentPatch(Centp);
  task->SetClusterTPC(ClusterTPC);
  task->SetDCAxy(DCAxy);
  task->SetFilterBit(FilterBit);
   //  Printf( "%f\n", Nsigma);
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("cHistLambda_%s",Dirsuffixname.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:lambdastar_%s", AliAnalysisManager::GetCommonFileName(),Dirsuffixname.Data()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("cHistPrim_%s",Dirsuffixname.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:primary_%s", AliAnalysisManager::GetCommonFileName(),Dirsuffixname.Data()));

	mgr->ConnectInput (task, 0, cinput0);
	mgr->ConnectOutput(task,1,coutput1);
	mgr->ConnectOutput(task,2,coutput2);

  return task;
}
