// #include "Config_amichalik_LMEE_MC.C"
// #include "LMEECutLib_amichalik.C"

AliAnalysisTask *AddTask_amichalik_MC2(Char_t* outputFileName="LMEEoutput.root",
                                              Bool_t getFromAlien=kFALSE,
                                              Int_t triggerNames=(AliVEvent::kINT7),
                                              Int_t collCands=AliVEvent::kINT7 )
{
  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTast_amichalik_MC2", "No analysis manager found.");
    return 0;
  }
  cout << "NEW CLASSES USED" << std::endl;
  //  create task and add it to the manager
  	// gSystem->AddIncludePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE");

  TString configBasePath("$TRAIN_ROOT/caklein_lowmass/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";


  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/Config_amichalik_LMEE_MC.C ."))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/LMEECutLib_amichalik.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFile("Config_amichalik_LMEE_MC.C");
  TString configLMEECutLib("LMEECutLib_amichalik.C");

  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  //AOD Usage currently tested with separate task, to be merged
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011AOD","switching on ESD specific code");
    bESDANA=kTRUE;
  }


  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  //load dielectron configuration files
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())){
    gROOT->LoadMacro(configLMEECutLibPath.Data());
    cout << "Cut library loaded successfully" << endl;
  }
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data())){
    gROOT->LoadMacro(configFilePath.Data());
    cout << "Config file loaded successfully" << endl;
  }

  cout << "Everything should be loaded by now" << endl;

  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerNames); // Should be not doing ANYTHING!!!!
  task->SelectCollisionCandidates(collCands);

  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kStandard));
  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!


  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_amichalik_LMEEPbPb(i,hasMC,bESDANA);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("amichalik_LMEEPbPb_tree",
                       TTree::Class(),
                       AliAnalysisManager::kExchangeContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("amichalik_LMEEPbPb_out",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("amichalik_LMEEPbPb_CF",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("amichalik_EventStatPbPb",
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
