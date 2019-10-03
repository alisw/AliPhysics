AliAnalysisTask *AddTask_Remi_LMEEPbPb2011AOD_SemiCent2(Char_t* outputFileName="LMEEoutput.root", Bool_t runAll=kFALSE,Bool_t setMC=kFALSE,Bool_t getFromAlien=kFALSE, Bool_t PIDbaseline=kFALSE, Bool_t rejOnly=kTRUE) {
  
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
  
  
  if (rejOnly) {
    
    //    AliDielectron *lowmass1=ConfigRemiLMEEPbPb2011AOD(1,hasMC,bESDANA);
    //    task->AddDielectron(lowmass1);
    //    printf("add: %s\n",lowmass1->GetName());
    
    
    //    AliDielectron *lowmass2=ConfigRemiLMEEPbPb2011AOD(2,hasMC,bESDANA);
    //    task->AddDielectron(lowmass2);
    //    printf("add: %s\n",lowmass2->GetName());
    
  
      AliDielectron *lowmass3=ConfigRemiLMEEPbPb2011AOD(3,hasMC,bESDANA);
      task->AddDielectron(lowmass3);
      printf("add: %s\n",lowmass3->GetName());
      AliDielectron *lowmass7=ConfigRemiLMEEPbPb2011AOD(7,hasMC,bESDANA);
      task->AddDielectron(lowmass7);
      printf("add: %s\n",lowmass7->GetName());

      AliDielectron *lowmass11=ConfigRemiLMEEPbPb2011AOD(11,hasMC,bESDANA);
      task->AddDielectron(lowmass11);
      printf("add: %s\n",lowmass11->GetName());

      AliDielectron *lowmass15=ConfigRemiLMEEPbPb2011AOD(15,hasMC,bESDANA);
      task->AddDielectron(lowmass15);
      printf("add: %s\n",lowmass15->GetName());

      AliDielectron *lowmass19=ConfigRemiLMEEPbPb2011AOD(19,hasMC,bESDANA);
      task->AddDielectron(lowmass19);
      printf("add: %s\n",lowmass19->GetName());
      AliDielectron *lowmass23=ConfigRemiLMEEPbPb2011AOD(23,hasMC,bESDANA);
      task->AddDielectron(lowmass23);
      printf("add: %s\n",lowmass23->GetName());
      AliDielectron *lowmass27=ConfigRemiLMEEPbPb2011AOD(27,hasMC,bESDANA);
      task->AddDielectron(lowmass27);
      printf("add: %s\n",lowmass27->GetName());


  
  //  AliDielectron *lowmass4=ConfigRemiLMEEPbPb2011AOD(4,hasMC,bESDANA);
  //  task->AddDielectron(lowmass4);
  //  printf("add: %s\n",lowmass4->GetName());

  
  //  AliDielectron *lowmass5=ConfigRemiLMEEPbPb2011AOD(5,hasMC,bESDANA);
  //  task->AddDielectron(lowmass5);
  //  printf("add: %s\n",lowmass5->GetName());


}
else {
  AliDielectron *lowmass4=ConfigRemiLMEEPbPb2011AOD(4,hasMC,bESDANA);
  task->AddDielectron(lowmass4);
  printf("add: %s\n",lowmass4->GetName());


  AliDielectron *lowmass1=ConfigRemiLMEEPbPb2011AOD(1,hasMC,bESDANA);
  task->AddDielectron(lowmass1);
  printf("add: %s\n",lowmass1->GetName());


if (PIDbaseline) {
	AliDielectron *lowmass7=ConfigRemiLMEEPbPb2011AOD(7,hasMC,bESDANA);
	task->AddDielectron(lowmass7);
	printf("add: %s\n",lowmass7->GetName());

 }
 }

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent2_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent2_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("Remi_LMEEPbPb2011SemiCent2_CF",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("Remi_EventStatPbPb2011SemiCent2",
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
