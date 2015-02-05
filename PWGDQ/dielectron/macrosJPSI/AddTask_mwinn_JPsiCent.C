AliAnalysisTask *AddTask_mwinn_JPsiCent(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_JPsi", "No analysis manager found.");
    return 0;
  }


  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  
  //set config file name
  //  TString configFile("$TRAIN_ROOT/mwinn_jpsiCorr/ConfigJpsi2eeDataCent_light.C");
  TString configFile("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mw_pPbCENT.C");
  //  TString configFile("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mw_pPbcent.C");
  //TString configFile("$TRAIN_ROOT/jpsi_JPSI/ConfigJpsiStandard.C");
  if ( trainConfig.Contains("PbPb") ) configFile="$TRAIN_ROOT/jpsi_JPSI/ConfigJpsi2eePbPb.C";
  
  //if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
  //  ::Info("AddTaskJPSI", "Using AOD configuration");
  //  configFile="$TRAIN_ROOT/util/dielectron/dielectron/macros/ConfigJpsi2eeDataAOD.C";
  // }

  TString list=gSystem->Getenv("LIST");
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDieData");
   if (!hasMC&&(!list.Contains("LHC12h")) ) task->UsePhysicsSelection();//taking out for testing//taken out for AOD
   task->SetTriggerMask(AliVEvent::kTRD|AliVEvent::kINT7);
  /*  if (list.Contains("LHC13b"))task->SetTriggerMask(AliVEvent::kINT7); //kINT7?, MB-trigger for MB pPb
  if (list.Contains("LHC13c"))task->SetTriggerMask(AliVEvent::kINT7); //kINT7?, MB-trigger for MB pPb
  if (list.Contains("LHC13d"))task->SetTriggerMask(AliVEvent::kTRD); //TRD trigger
  if (list.Contains("LHC13e"))task->SetTriggerMask(AliVEvent::kTRD); //TRD trigger//only for testing
  if (list.Contains("LHC11d")) task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);
  //if (list.Contains("LHC12h")) task->SetTRDtrigger(1+2);*/
  mgr->AddTask(task);

  
  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  cout << "!!!!!!!!!!!! nDie = " << nDie << endl;
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi2ee(i);
    if (!jpsi) continue;
    task->AddDielectron(jpsi);
//    printf("add: %s\n",jpsi->GetName());
  }
  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  //********************************
  //AOD selection not working, because default Vertex in Contructor of AliDielectronEventCuts not available for AODs!!!!!!!
  //Therefore chose default one for AODs
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  //******************
  eventCuts->SetRequireVertex();//NOTE: all of these cuts can for some reasons not be applied to self-filtered AODs by mwinn in, also true for 2013AODs bug in Dielectron-framework
  //(list/hera/alice/mwinn/mwinn/train/lists/...)
  //  eventCuts->SetRequire2013vertexandevent();//taken out for AOD
  eventCuts->SetMinVtxContributors(1);//taken out for AOD
  eventCuts->SetVertexZ(-10.,10.);//taken out for AOD
  //  eventCuts->Print();
//   eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTracksOrSPD);
//   eventCuts->SetRequireV0and();
//   if ( trainConfig=="PbPb" ){
//     eventCuts->SetCutOnMultipicityITSTPC();
//   }
   task->SetEventFilter(eventCuts);

//   task->SetTriggerOnV0AND();
  if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("jpsi_mwinn_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "jpsi_mwinn_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsiAOD_mwinn_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_mwinn_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_mwinn_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
