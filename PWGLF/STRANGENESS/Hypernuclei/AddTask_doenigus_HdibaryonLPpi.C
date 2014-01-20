AliAnalysisTask *AddTask_doenigus_HdibaryonLPpi(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_doenigus_HdibaryonLPpi", "No analysis manager found.");
    return 0;
  }
  /*
  //============= Set Task Name ===================
  TString taskName=("AliAnalysisTaskHdibaryonLPpi.cxx+");
  //===============================================
  //            Load the task
    
  if (gProof){
    TString taskSO=gSystem->pwd();
    taskSO+="/";
    taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
    gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
  }
  
  gROOT->LoadMacro("AliAnalysisTaskHdibaryonLPpi.cxx+");
  gROOT->LoadMacro(taskName.Data());
  */

   Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

   // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   //Bool_t isMC=kFALSE; // kTRUE in case of MC
   //if (hasMC) isMC=kTRUE;
   //AddTaskPIDResponse(isMC); 

   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   Bool_t isMC=kFALSE; // kTRUE in case of MC
   if (hasMC) isMC=kTRUE;
   Bool_t useTPCEtaCorrection = kFALSE;
   AddTaskPIDResponse(isMC, kTRUE, kFALSE, 2, kFALSE, "", useTPCEtaCorrection);

  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskSE *taskHdibaryonLPpi = new AliAnalysisTaskHdibaryonLPpi("doenigus_HdibaryonLPpi");
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  // taskHdibaryonLPpi->SelectCollisionCandidates(AliVEvent::kAny);
  //  taskTest->SelectCollisionCandidates(AliVEvent::kMB);

  mgr->AddTask(taskHdibaryonLPpi);

  if (!hasMC){
  //  taskLambda1520MC->SelectCollisionCandidates(AliVEvent::kMB);
  taskHdibaryonLPpi->SelectCollisionCandidates(AliVEvent::kAny);
  }
  //  taskTest->SelectCollisionCandidates(AliVEvent::kMB);
  //  mgr->AddTask(taskLambda1520MC);

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //this is the old way!!!
  //AliAnalysisDataContainer *cinput  = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");

  //AliAnalysisDataContainer *coutput0_doenigusHdibaryonLPpi = mgr->CreateContainer("doenigus_tree_doenigusHdibaryonLPpi", TTree::Class(), AliAnalysisManager::kExchangeContainer, "doenigus_default_doenigusHdibaryonLPpi");  

  //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput_doenigusHdibaryonLPpi = mgr->CreateContainer("doenigus_HdibaryonLPpi", TList::Class(), AliAnalysisManager::kOutputContainer,"doenigus_HdibaryonLPpi.root");

  //           connect containers
  mgr->ConnectInput  (taskHdibaryonLPpi,  0, cinput );
  //mgr->ConnectOutput (taskHdibaryonLPpi,  0, coutput0_doenigusHdibaryonLPpi); 
  mgr->ConnectOutput (taskHdibaryonLPpi,  1, coutput_doenigusHdibaryonLPpi);

  return taskHdibaryonLPpi;
}
