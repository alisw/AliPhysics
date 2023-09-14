void AddTask_GammaPureMC(Int_t isK0 = 0 , Double_t maxpT = 100 , Int_t doMultStudies = 0, int doJetStudies = 0) {

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_GammaPureMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaPureMC *task=NULL;
  task= new AliAnalysisTaskGammaPureMC("GammaPureMC");
  // if no k0 desired, set to isK0
  task->SetIsK0(isK0);
  task->SetMaxPt(maxpT);
  task->SetDoMultStudies(doMultStudies);
  task->SetDoJetStudies(doJetStudies);

  TString AddString = "";
  if(doMultStudies){
    AddString+=Form("_Mult%i", doMultStudies);
  }
  if(doJetStudies){
    AddString+=Form("_Jet%i", doJetStudies);
  }
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaPureMC%s", AddString.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaPureMC%s",AliAnalysisManager::GetCommonFileName(), AddString.Data()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
