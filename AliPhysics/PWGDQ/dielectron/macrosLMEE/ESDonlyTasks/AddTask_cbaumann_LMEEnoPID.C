AliAnalysisTask *AddTask_cbaumann_LMEEnoPID(Bool_t withMC = kFALSE,Bool_t enableCF = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
	Error("AddTask_cbaumann_LMEEnoPID", "No analysis manager found.");
	return 0;
  }

  //create config File names: TRAIN_ROOT is for running on GSI train, 
  // ALICE_ROOT for CERN Lego trains
  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";
  TString configFile("ConfigLMEEnoPID.C");
  TString configLMEECutLib("LMEECutLib.C");

  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);


  //AOD Usage not yet testes/avialable-------------------------------------

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
	::Info("AddTaskLMEEnoPID", "no dedicated AOD configuration");
//	configFile="$TRAIN_ROOT/util/dielectron/dielectron/macros/ConfigLMEEnoPIDAOD.C";	
	
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if (!hasMC) hasMC=withMC;

//  create task and add it to the manager

  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
	gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
	gROOT->LoadMacro(configFilePath.Data());


  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kpp2010TPCandTOF)); //
	

  //load dielectron configuration file

  //add dielectron analysis with different cuts to the task
   //no memleak
  AliDielectron *lowmass0=ConfigLMEEnoPID(0,hasMC,enableCF);
  task->AddDielectron(lowmass0);
  printf("add: %s\n",lowmass0->GetName());

  AliDielectron *lowmass1=ConfigLMEEnoPID(1,hasMC,enableCF);
  task->AddDielectron(lowmass1);
  printf("add: %s\n",lowmass1->GetName());

   AliDielectron *lowmass2=ConfigLMEEnoPID(2,hasMC,enableCF);
  task->AddDielectron(lowmass2);
  printf("add: %s\n",lowmass2->GetName());
   
  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("cbaumann_LMEEnoPID_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		"cbaumann_LMEEnoPID_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("cbaumann_LMEEnoPID_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEnoPID_out.root");

  AliAnalysisDataContainer *cOutputHist2 = 0x0;
 if (enableCF) {
	cOutputHist2 =
	mgr->CreateContainer("cbaumann_LMEEnoPID_CF",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEnoPID_out.root");
 }
  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("cbaumann_EventStatnoPID",
		TH1D::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEnoPID_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
if (enableCF)  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
