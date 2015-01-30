AliAnalysisTask *AddTask_sweber_JPsi_pPb_TRDtrigger(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sweber_JPsi_pPb_TRDtrigger", "No analysis manager found.");
    return 0;
  }


  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  
  //set config file name
  // TString configFile("$TRAIN_ROOT/mwinn_jpsiCorr/ConfigJpsi2eeData.C");
  TString configFile("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_sweber_pPb_TRDtrigger.C");
  //TString configFile("$TRAIN_ROOT/jpsi_JPSI/ConfigJpsiStandard.C");
  // if ( trainConfig.Contains("PbPb") ) configFile="$TRAIN_ROOT/jpsi_JPSI/ConfigJpsi2eePbPb.C";
  
  //if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
  //  ::Info("AddTaskJPSI", "Using AOD configuration");
  //  configFile="$TRAIN_ROOT/util/dielectron/dielectron/macros/ConfigJpsi2eeDataAOD.C";
  // }

  TString list=gSystem->Getenv("LIST");
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDieDataTRDtrigger");
  if (!hasMC ) task->UsePhysicsSelection();//taking out for testing
  task->SetTriggerMask(AliVEvent::kTRD); 
  //task->SetTRDtrigger(3);
  task->SetFiredTriggerName("CINT7WUHJT-B-NOPF-CENT",kTRUE);//take not kTRUE in order not to get both triggers, exclude (since not yet otherwise possible)
  //
  //not yet implemented to get a logical or of 2 different trigger classes, therefore exclude HJT, not exactly what I want...	
  //   if (list.Contains("LHC11d")) task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);
  //if (list.Contains("LHC12h")) task->SetTRDtrigger(1+2);
  mgr->AddTask(task);

  
  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  cout << "!!!!!!!!!!!! nDie = " << nDie << endl;
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_sweber_pPb_TRDtrigger(i);
    if (!jpsi) continue;
    task->AddDielectron(jpsi);
//    printf("add: %s\n",jpsi->GetName());
  }
  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();//NOTE: all of these cuts can for some reasons not be applied to self-filtered AODs by mwinn in 
  //(list/hera/alice/mwinn/mwinn/train/lists/...)
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
//   eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTracksOrSPD);
//   eventCuts->SetRequireV0and();
//   if ( trainConfig=="PbPb" ){
//     eventCuts->SetCutOnMultipicityITSTPC();
//   }
   task->SetEventFilter(eventCuts);

//   task->SetTriggerOnV0AND();
  if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container
  TString containerName = "AnalysisResults.root";
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("sweber_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("sweber_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("sweber_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("sweber_EventStat",
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
