AliAnalysisTask *AddTask_miweber_LMEE_pp(TString outputFileName = "AnalysisResult.root", TString directoryBaseName = "miweber_LMEE_pp"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_miweber_LMEE_pp", "No analysis manager found.");
    return 0;
  }

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler

  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_miweber_LMEE_pp.C");
  TString configLMEECutLib("LMEECutLib_miweber.C");
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);
  
   
  //load dielectron configuration files
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());

  // cut lib
  LMEECutLib* cutlib = new LMEECutLib();

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AOD Usage currently tested with Input handler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_pp", "no dedicated AOD configuration");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_pp","switching on ESD specific code");
    bESDANA=kTRUE;
  }
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData_miweber_pp");
  if (!hasMC) task->UsePhysicsSelection();

  //Add event filter
  Int_t triggerNames=(AliVEvent::kINT7+AliVEvent::kMB+AliVEvent::kINT8);//pp 2010/2011 Min Bias?
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kpp2010All));
  task->SelectCollisionCandidates(triggerNames);
  task->SetTriggerMask(triggerNames);
  task->SetRejectPileup();
  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!

  // Add the task to the manager
  mgr->AddTask(task);
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_miweber_pp(i,hasMC,bESDANA);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop


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
  //                         "miweber_LMEE_pp_CF.root");
  
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
