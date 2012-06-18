AliAnalysisTask *AddTask_mwinn_JPsi(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_mwinn_JPsi", "No analysis manager found.");
    return 0;
  }

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTask_mwinn_JPsi", "Not used in AOD!!!");
    return 0;
  }

  //Get the current train configuration
  //NOTE train config not there on grid, PbPb anyhow not used for correlations
  //TString trainConfig=gSystem->Getenv("CONFIG_FILE");

  //set config file name
  //for grid running: "$ALICE_ROOT/PWGDQ/dielectron/macros/ConfigJpsi2eeDebugTree.C"
  //for gsi running "$TRAIN_ROOT/jpsi_JPSI/ConfigJpsi2eeDebugTree.C"
  TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mw_pp.C");
  //NOTE train config not there on grid, PbPb anyhow not used for correlations
  //if ( trainConfig.Contains("PbPb") ) configFile="$TRAIN_ROOT/jpsi_JPSI/ConfigJpsi2eeDebugTreePbPb.C";
  
    //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());

  
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("mwinn_JPsi");
  if (!hasMC) task->UsePhysicsSelection();
	  //NOTE trigger selection for EMCAL triggered data...
	task->SetTriggerMask(AliVEvent::kEMC7);
  mgr->AddTask(task);
  
  //load dielectron configuration file, already done before!
  //gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
	//NOTE: for 11d running just taking one single configuration, see configjpsi2eeDebugTree
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_mw_pp(i);
    if(jpsi) {
      jpsi->SetDontClearArrays();
      task->AddDielectron(jpsi);
    }
  }

  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //NOTE train config not there on grid, PbPb anyhow not used for correlations
  //  if ( trainConfig.Contains("PbPb") ) eventCuts->SetCentralityRange(40.,80.);
  task->SetEventFilter(eventCuts);
  
  //   task->SetTriggerOnV0AND();
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("mwinn_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "mwinn_default");
  												 
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("mwinn_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "mwinn.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("mwinn_CF",
                         TList::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "mwinn.root");
												 


  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("mwinn_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "mwinn.root");

												 
	//is the first output type  a problem??!
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  
  return task;
}
