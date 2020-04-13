AliAnalysisTask *AddTask_LNNv0Bkg(Int_t bkg, Bool_t isMC){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_LNNv0Bkg", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====
 TString name(Form("LNN_Bkg_Type_%i",bkg));

  AliAnalysisTaskLNNv0Bkg *taskLNNv0Bkg = new AliAnalysisTaskLNNv0Bkg(name.Data());
  taskLNNv0Bkg->SetBkgType(bkg);
  if(isMC) taskLNNv0Bkg->SetMC(); 
  mgr->AddTask(taskLNNv0Bkg);
  
  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  
  TString outputFileName("LNNv0Bkg.root"); /*AliAnalysisManager::GetCommonFileName();*/

 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("LNNlist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
   
  //           connect containers
  mgr->ConnectInput  (taskLNNv0Bkg,  0, cinput );
  mgr->ConnectOutput (taskLNNv0Bkg,  1, coutput1);

  return taskLNNv0Bkg;
}
