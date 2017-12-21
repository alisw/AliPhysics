AliAnalysisTask *AddTask_miweber_LMEE_V0(Int_t cutDefinition = 0, TString outputFileName = "AnalysisResult.root", TString directoryBaseName = "miweber_LMEE_V0", Bool_t isNano = kFALSE, Bool_t bCutQA = kTRUE, Bool_t bCutV0 = kTRUE, Bool_t bCutV0Tight = kFALSE){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_miweber_LMEE_V0", "No analysis manager found.");
    return 0;
  }

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler

  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_miweber_LMEE_V0.C");
  TString configFilePath(configBasePath+configFile);
  
   
  //load dielectron configuration files
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AOD Usage currently tested with Input handler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_V0", "no dedicated AOD configuration");
    return NULL;
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_V0","switching on ESD specific code");
    bESDANA=kTRUE;
  }
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiElData_miweber_V0");

  // for MC no  need for physics selection and for Nano AODs this has been done already  
  if (!hasMC && !isNano) task->UsePhysicsSelection();

  //Add event filter
  Int_t triggerNames = AliVEvent::kINT7;//PbPb Min Bias, can be set also from outside

  // for Nano AODs this has been done already  
  if(!isNano){
    task->SelectCollisionCandidates(triggerNames);
    task->SetTriggerMask(triggerNames);
    // task->SetRejectPileup(); // to be done differently (too strong cuts at the moment in dielectron framework) 
  }

  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!

  //Add event filter
  task->SetEventFilter( GetEventCuts() );

  // Add the task to the manager
  mgr->AddTask(task);
  
  //add dielectron analysis with selected cut to the task
  AliDielectron *diel_low = Config_miweber_LMEE_V0(cutDefinition,bESDANA,bCutQA,bCutV0,bCutV0Tight);
  if(diel_low){

    // NOT NEEDED HERE:
    // AliDielectronVarCuts *eventplaneCuts = new AliDielectronVarCuts("eventplaneCuts","eventplaneCuts");
    // eventplaneCuts->AddCut(AliDielectronVarManager::kQnTPCrpH2,-999.,kTRUE); // makes sure that the event has an eventplane
    // eventplaneCuts->Print();
    // diel_low->GetEventFilter().AddCuts(eventplaneCuts);
  
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("%s_tree",directoryBaseName.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("%s_out",directoryBaseName.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("%s_CF",directoryBaseName.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  //                         "miweber_LMEE_V0_CF.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("%s_EventStat",directoryBaseName.Data()),
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
