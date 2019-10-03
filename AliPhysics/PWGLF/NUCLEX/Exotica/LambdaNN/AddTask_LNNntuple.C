AliAnalysisTask *AddTask_LNNntuple(Bool_t isMC=kFALSE){

 //get the current analysis manager
 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
  Error("AddTask_LH3Pi", "No analysis manager found.");
  return 0;
 }

 // mc event handler                                                                                                                                                                                                           
 AliMCEventHandler* mchandler = new AliMCEventHandler();
 // Not reading track references                                                                                                                                                                                           
 mchandler->SetReadTR(kFALSE);
 mgr->SetMCtruthEventHandler(mchandler);

 //========= Add task to the ANALYSIS manager =====

 AliAnalysisTaskLNNntuple *taskLNN = new AliAnalysisTaskLNNntuple("LNNntuple",isMC);

 if(!taskLNN) printf("\n\n\n   *************** NO TASK LOADED!! ************ \n\n\n\n");

 mgr->AddTask(taskLNN);

 //================================================
 //              data containers
 //================================================
 //            find input container
 TString outputFileName = AliAnalysisManager::GetCommonFileName();
 printf("Output filename %s \n",outputFileName.Data());
 AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
 AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("LNNlist", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());  

 //           connect containers
 mgr->ConnectInput  (taskLNN,  0, cinput );
 mgr->ConnectOutput (taskLNN,  1, coutput1);

 return taskLNN;
}
