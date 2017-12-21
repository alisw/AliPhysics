AliAnalysisTask *AddTask_caklein_LMEEPbPb_AOD(Char_t* outputFileName="LMEEoutput.root",
                                              TString configFile="Config_caklein_LMEEPbPb_AOD.C",
                                              Bool_t getFromAlien=kFALSE,
                                              Bool_t cutlibPreloaded=kFALSE,
                                              Int_t triggerNames=AliVEvent::kINT7,
                                              Int_t collCands=AliVEvent::kINT7,
                                              Int_t eventCut=0,
                                              Int_t wagonnr=0,
                                              Int_t centrality=4 /* 3=minbias, 2=50-80, 1=10-50, 0=00-10 */)
{
  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_caklein_LMEEPbPb_AOD", "No analysis manager found.");
    return 0;
  }
  cout << "NEW CLASSES USED" << std::endl;
  //  create task and add it to the manager
  //	gSystem->AddIncludePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE");

  TString configBasePath("$TRAIN_ROOT/caklein_lowmass/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";


  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/%s .", configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/LMEECutLib_caklein.C ."))
      ) {
    cout << "Copy config from Alien" << std::endl;
    configBasePath=Form("%s/",gSystem->pwd());
  }

  // TString configFile("Config_caklein_LMEEPbPb_AOD.C");
  TString configLMEECutLib("LMEECutLib_caklein.C");
  TString configFilePath(configBasePath+configFile.Data());
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  Bool_t err=kFALSE;
  if (!cutlibPreloaded) { // should not be needed but seems to be...
    std::cout << "Cutlib was not preloaded" << std::endl;
    err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
    err |= gROOT->LoadMacro(configFilePath.Data());
  }
  else{
    std::cout << "Cutlib was preloaded in a previous task" << std::endl;
  }
  // err |= gROOT->LoadMacro(configFilePath.Data());
  if (err) { Error("AddTask_caklein_LMEEPbPb_AOD","Config(s) could not be loaded!"); return 0x0; }

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011AOD","switching on ESD specific code");
    bESDANA=kTRUE;
  }


  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  cout << "Everything should be loaded by now" << endl;

  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron(Form("MultiDiEData%d", wagonnr));
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerNames); // Should be not doing ANYTHING!!!!
  task->SelectCollisionCandidates(collCands);

  task->SetEventFilter(cutlib->GetEventCuts(eventCut));
  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!


  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_caklein_LMEEPbPb_AOD(arrNames->At(i)->GetName(),hasMC,bESDANA,centrality);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer(Form("caklein_LMEEPbPb_tree%d", wagonnr),
                       TTree::Class(),
                       AliAnalysisManager::kExchangeContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer(Form("caklein_LMEEPbPb_out%d", wagonnr),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer(Form("caklein_LMEEPbPb_CF%d", wagonnr),
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer(Form("caklein_EventStatPbPb%d", wagonnr),
                       TH1D::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
