AliAnalysisTask *AddTask_caklein_LMEEPbPb_AOD_2018(Char_t* outputFileName = "LMEE.root",
                                              TString configFile = "Config_caklein_LMEEPbPb_AOD_2018.C",
                                              TString configLMEECutLib = "LMEECutLib_caklein_2018.C",
                                              Bool_t getFromAlien    = kFALSE, // load config from AliEn
                                              Int_t collCands = AliVEvent::kINT7,
                                              Int_t eventCut = 0, // use different event cuts
                                              Bool_t cutlibPreloaded = kFALSE,
                                              Int_t wagonnr  = 0, // needs to be set when using different wagons per trains
                                              Int_t centrality = 4 /* 3=minbias, 2=50-80, 1=10-50, 0=00-10 */)
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_caklein_LMEEPbPb_AOD", "No analysis manager found.");
    return 0;
  }

  // If used at GSI use code lying there
  TString configBasePath("$TRAIN_ROOT/caklein_lowmass/");
  TString trainRoot = gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";


  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/%s .", configFile.Data())))
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/%s .", configLMEECutLib.Data())))
      ) {
    cout << "Copy config from Alien" << std::endl;
    configBasePath=Form("%s/",gSystem->pwd());
  }
  
  // TString configFile("Config_caklein_LMEEPbPb_AOD.C");
  TString configFilePath      (configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  // Load cut library when running several trains at once. Library has to be only loaded once otherwise the trains fail miserably
  Bool_t err=kFALSE;
  if (!cutlibPreloaded) { // should not be needed but seems to be...
    std::cout << "Cutlib was not preloaded and will be loaded now" << std::endl;
    err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
    err |= gROOT->LoadMacro(configFilePath.Data());
  }
  else{
    std::cout << "Cutlib was preloaded in a previous task, no need to load it again" << std::endl;
  }
  if (err) { Error("AddTask_caklein_LMEEPbPb_AOD","Config(s) could not be loaded!"); return 0x0; }

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  std::cout << "Everything should be loaded by now" << std::endl;

  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron(Form("MultiDiEData%d", wagonnr));
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(collCands); // Should be not doing ANYTHING!!!! Introduced by Mahmut
  task->SelectCollisionCandidates(collCands);

  task->SetEventFilter(cutlib->GetEventCuts(eventCut));

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i < nDie; ++i){ //nDie defined in config file
    AliDielectron *diel_low = Config_caklein_LMEEPbPb_AOD_2018(arrNames->At(i)->GetName(),hasMC,centrality);
    if(!diel_low) continue;
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
  mgr->ConnectOutput(task, 1, cOutputHist1);

  return task;
}
