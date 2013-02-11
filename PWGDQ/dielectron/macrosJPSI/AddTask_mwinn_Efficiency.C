AliAnalysisTask *AddTask_mwinn_Efficiency(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_Efficiency", "No analysis manager found.");
    return 0;
  }

  //set config file name
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTask_jpsi_Efficiency", "Not running in AOD");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDieDataEfficiency");
  if (!hasMC) task->UsePhysicsSelection();
  mgr->AddTask(task);
  
  //load dielectron configuration file
  TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mw_EffpPb.C");
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDieEff; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_mw_EffpPb(i);
    task->AddDielectron(jpsi);
  }
  

  //create output container
  TString containerName = "JPSI.root";
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("mwinnEff_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
			 containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("mwinnEff_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("mwinnEff_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
			 containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("mwinnEff_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
			 containerName.Data());
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
