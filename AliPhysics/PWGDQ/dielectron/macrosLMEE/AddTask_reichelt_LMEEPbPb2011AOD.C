AliAnalysisTask *AddTask_reichelt_LMEEPbPb2011AOD(Bool_t getFromAlien=kFALSE,
                                                  TString configFile="Config_reichelt_LMEEPbPb2011.C",
                                                  Bool_t cutlibPreloaded=kFALSE, Bool_t flag1=kFALSE,
                                                  Char_t* outputFileName="LMEEoutput.root",
                                                  Int_t triggerNames=(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral),
                                                  Int_t collCands=AliVEvent::kAny)
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLMEEPbPb2011AOD", "No analysis manager found.");
    return 0;
  }
  
  // environment for testing/running at GSI:
  TString configBasePath("$TRAIN_ROOT/reichelt_lowmass/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  // typical Aliroot environment:
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  ::Info("AddTaskLMEEPbPb2011AOD",Form("configBasePath.Data(): %s\n",configBasePath.Data()));
  
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
      && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/p/preichel/PWGDQ/dielectron/macrosLMEE/LMEECutLib_reichelt.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  
  TString configLMEECutLib("LMEECutLib_reichelt.C");
  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);
  
  //load dielectron configuration files
  Bool_t err=kFALSE;
  if (!cutlibPreloaded) { // should not be needed but seems to be...
    //if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())) ///...because this check doesnt work
    err |= gROOT->LoadMacro(configLMEECutLibPath.Data());
  }
  err |= gROOT->LoadMacro(configFilePath.Data());
  if (err) { Error("AddTaskLMEEPbPb2011AOD","Config(s) could not be loaded!"); return 0x0; }
  
  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011AOD","running on AODs.");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011AOD","switching on ESD specific code, make sure ESD cuts are used.");
    bESDANA=kTRUE;
  }
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  std::cout << "hasMC = " << hasMC << std::endl;
  
  // Set up the task
  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  if (!hasMC) task->SetTriggerMask(triggerNames);
  task->SelectCollisionCandidates(collCands);  
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1, hasMC));
  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!
  task->SetRandomizeDaughters(randomizeDau);//default kFALSE
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *diel_low = Config_reichelt_LMEEPbPb2011(i, hasMC, bESDANA);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop
  
  mgr->AddTask(task);
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("reichelt_LMEEPbPb2011_tree",
                       TTree::Class(),
                       AliAnalysisManager::kExchangeContainer,
                       outputFileName);
  
  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("reichelt_LMEEPbPb2011_out",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);
  
  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("reichelt_LMEEPbPb2011_CF",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);
  
  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("reichelt_EventStatPbPb2011",
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
