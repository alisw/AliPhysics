AliAnalysisTask *AddTask_Remi_LMEEPbPb2011AOD_SemiCent1(Char_t* outputFileName="LMEEoutput.root", Bool_t runAll=kFALSE,Bool_t setMC=kFALSE,Bool_t getFromAlien=kFALSE, Bool_t PIDbaseline=kFALSE, Bool_t rejOnly=kTRUE) {
  
  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_Remi_LMEEPbPb2011", "No analysis manager found.");
    return 0;
  }
  
  
  //  create task and add it to the manager
  //	gSystem->AddIncludePath("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE");
  
  
  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  //  if (trainRoot.IsNull()) configBasePath= "/home/tanizaki/nfs/LMee_Deflection/ver1/";
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";
  
  
  if (getFromAlien &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/ConfigRemiLMEEPbPb2011AOD.C .")) &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/LMEECutLibRemi.C ."))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }
  
  
  TString configFile("ConfigRemiLMEEPbPb2011AOD.C");
  TString configLMEECutLib("LMEECutLibRemi.C");
  
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
  Bool_t hasMC=setMC;
  if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0)
    hasMC=kTRUE;
  
  
  
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());
  
  
  LMEECutLibRemi* cutlib = new LMEECutLibRemi();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  //  task->SelectCollisionCandidates(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  //  task->SetRejectPileup();
  task->SelectCollisionCandidates(AliVEvent::kAny);  
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLibRemi::kPbPb2011pidITSTPCTOF)); //
  
  //load dielectron configuration file

  //add dielectron analysis with different cuts to the task

  AliDielectron *lowmass2=ConfigRemiLMEEPbPb2011AOD(2,hasMC,bESDANA);
  task->AddDielectron(lowmass2);
  printf("add: %s\n",lowmass2->GetName());
  AliDielectron *lowmass6=ConfigRemiLMEEPbPb2011AOD(6,hasMC,bESDANA);
  task->AddDielectron(lowmass6);
  printf("add: %s\n",lowmass6->GetName());
  AliDielectron *lowmass10=ConfigRemiLMEEPbPb2011AOD(10,hasMC,bESDANA);
  task->AddDielectron(lowmass10);
  printf("add: %s\n",lowmass10->GetName());
  AliDielectron *lowmass14=ConfigRemiLMEEPbPb2011AOD(14,hasMC,bESDANA);
  task->AddDielectron(lowmass14);
  printf("add: %s\n",lowmass14->GetName());

  AliDielectron *lowmass18=ConfigRemiLMEEPbPb2011AOD(18,hasMC,bESDANA);
  task->AddDielectron(lowmass18);
  printf("add: %s\n",lowmass18->GetName());
  AliDielectron *lowmass22=ConfigRemiLMEEPbPb2011AOD(22,hasMC,bESDANA);
  task->AddDielectron(lowmass22);
  printf("add: %s\n",lowmass22->GetName());
  AliDielectron *lowmass26=ConfigRemiLMEEPbPb2011AOD(26,hasMC,bESDANA);
  task->AddDielectron(lowmass26);
  printf("add: %s\n",lowmass26->GetName());


  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent1_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent1_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent1_CF",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("Remi_EventStatPbPb2011SemiCent1",
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
