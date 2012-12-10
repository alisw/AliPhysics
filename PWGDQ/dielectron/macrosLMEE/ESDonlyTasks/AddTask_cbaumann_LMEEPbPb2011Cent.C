AliAnalysisTask *AddTask_cbaumann_LMEEPbPb2011Cent(Bool_t runRejection=kFALSE, Bool_t setMC=kFALSE,Bool_t enableCF=kFALSE, Bool_t switchToPhiV=kTRUE,Bool_t switchToOA=kFALSE, Bool_t getFromAlien=kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
	Error("AddTask_cbaumann_LMEEPbPb2011", "No analysis manager found.");
	return 0;
  }

  // 3 options:
  // TRAIN_ROOT is for running on GSI train, 
  // ALICE_ROOT for CERN Lego trains
  // getFromAlien: load files from user ALIEN dir 
  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");                                                                            
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien && 
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/ConfigLMEEPbPb2011.C .")) &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/LMEECutLib.C ."))
     ) {
	configBasePath=Form("%s/",gSystem->pwd());
  }
  TString configFile("ConfigLMEEPbPb2011.C");
  TString configLMEECutLib("LMEECutLib.C");

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
  cutlib->SetMCFlag(hasMC);
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEDataPbPb2011Cent");
  if (!hasMC){ task->UsePhysicsSelection();
  }
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kPbPb2011TPCandTOF)); //

  //load dielectron configuration file

  //add dielectron analysis with different cuts to the task
  if (runRejection) {

	AliDielectron *lowmass4=ConfigLMEEPbPb2011(4,hasMC,enableCF);
	task->AddDielectron(lowmass4);
	printf("add: %s\n",lowmass4->GetName());
  }
  else {
	if (switchToPhiV) {

   AliDielectron *lowmass9=ConfigLMEEPbPb2011(9,hasMC,enableCF);
	  lowmass9->SetUseKF(kFALSE);
	  task->AddDielectron(lowmass9);
	  printf("add: %s\n",lowmass9->GetName()); }

	else if (switchToOA) {

   AliDielectron *lowmass12=ConfigLMEEPbPb2011(12,hasMC,enableCF);
	  lowmass12->SetUseKF(kFALSE);
	  task->AddDielectron(lowmass12);
	  printf("add: %s\n",lowmass12->GetName()); }

	else {
	  AliDielectron *lowmass1=ConfigLMEEPbPb2011(1,hasMC,enableCF);
	  lowmass1->SetUseKF(kFALSE);
	  task->AddDielectron(lowmass1);
	  printf("add: %s\n",lowmass1->GetName());
	}
	/*      AliDielectron *lowmass7=ConfigLMEEPbPb2011(7,hasMC,enableCF);
			task->AddDielectron(lowmass7);
			printf("add: %s\n",lowmass7->GetName());*/
  }

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("cbaumann_LMEEPbPb2011Cent_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		"cbaumann_LMEEPbPb2011Cent_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("cbaumann_LMEEPbPb2011Cent_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEPbPb2011Cent_out.root");

  AliAnalysisDataContainer *cOutputHist2 = 0x0;
  if (enableCF) {
	cOutputHist2 = 
	  mgr->CreateContainer("cbaumann_LMEEPbPb2011Cent_CF",
		  TList::Class(),
		  AliAnalysisManager::kOutputContainer,
		  "cbaumann_LMEEPbPb2011Cent_out.root");

  }
  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("cbaumann_EventStatPbPb2011",
		TH1D::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEPbPb2011Cent_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  if (enableCF) {
	mgr->ConnectOutput(task, 2, cOutputHist2);
  }
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
