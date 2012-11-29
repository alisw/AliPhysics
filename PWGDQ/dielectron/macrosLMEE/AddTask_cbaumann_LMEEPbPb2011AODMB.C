AliAnalysisTask *AddTask_cbaumann_LMEEPbPb2011AODMB(Bool_t runAll=kFALSE,Bool_t setMC=kFALSE,Bool_t getFromAlien=kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
	Error("AddTask_cbaumann_LMEEPbPb2011", "No analysis manager found.");
	return 0;
  }


//  create task and add it to the manager
//	gSystem->AddIncludePath("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE");


  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";


  if (getFromAlien &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/ConfigLMEEPbPb2011AOD.C .")) &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/LMEECutLibAOD.C ."))
     ) {
        configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFile("ConfigLMEEPbPb2011AOD.C");
  TString configLMEECutLib("LMEECutLibAOD.C");

  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  //AOD Usage currently tested with separate task, to be merged
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
	::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }

  //Do we have an MC handler?
  Bool_t hasMC=setMC;
  if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0)
	hasMC=kTRUE;



  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
	gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
	gROOT->LoadMacro(configFilePath.Data());


  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(AliVEvent::kMB);
  //task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  
task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kPbPb2011TPCandTOF)); //
	

  //load dielectron configuration file

  //add dielectron analysis with different cuts to the task
  AliDielectron *lowmass4=ConfigLMEEPbPb2011AOD(4,hasMC);
  task->AddDielectron(lowmass4);
  printf("add: %s\n",lowmass4->GetName());
/*
  AliDielectron *lowmass6=ConfigLMEEPbPb2011AOD(6,hasMC);
  task->AddDielectron(lowmass6);
  printf("add: %s\n",lowmass6->GetName());
*/
  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("cbaumann_LMEEPbPb2011_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		"cbaumann_LMEEPbPb2011_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("cbaumann_LMEEPbPb2011_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEPbPb2011_out.root");

  /*  AliAnalysisDataContainer *cOutputHist2 =
	  mgr->CreateContainer("cbaumann_lowmass_CF",
	  TList::Class(),
	  AliAnalysisManager::kOutputContainer,
	  "cbaumann_lowmass_CF.root");
	  */
  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("cbaumann_LMEEPbPb2011_CF",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEPbPb2011_out.root");

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("cbaumann_EventStatPbPb2011",
		TH1D::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEPbPb2011_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
