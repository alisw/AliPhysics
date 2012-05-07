AliAnalysisTask *AddTask_cbaumann_LMEEpp2010(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
	Error("AddTask_cbaumann_LMEEpp2010", "No analysis manager found.");
	return 0;
  }

  //set config file name
  TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/ConfigLMEEpp2010.C");

  //AOD Usage not yet testes/avialable-------------------------------------

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
	::Info("AddTaskLMEEpp2010", "no dedicated AOD configuration");
//	configFile="$TRAIN_ROOT/util/dielectron/dielectron/macros/ConfigLMEEpp2010AOD.C";	
	
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);


//  create task and add it to the manager

   gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/LMEECutLib.C");
  gROOT->LoadMacro(configFile.Data());
  LMEECutLib* cutlib = new LMEECutLib();
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kpp2010TPCandTOF)); //
	

  //load dielectron configuration file

  //add dielectron analysis with different cuts to the task
/*
  AliDielectron *lowmass1=ConfigLMEEpp2010(1,hasMC);
  task->AddDielectron(lowmass1);
  printf("add: %s\n",lowmass1->GetName());
*/
  AliDielectron *lowmass2=ConfigLMEEpp2010(2,hasMC);
  task->AddDielectron(lowmass2);
  printf("add: %s\n",lowmass2->GetName());


  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("cbaumann_LMEEpp2010_tree",
		TTree::Class(),
		AliAnalysisManager::kExchangeContainer,
		"cbaumann_LMEEpp2010_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("cbaumann_LMEEpp2010_out",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEpp2010_out.root");

  /*  AliAnalysisDataContainer *cOutputHist2 =
	  mgr->CreateContainer("cbaumann_lowmass_CF",
	  TList::Class(),
	  AliAnalysisManager::kOutputContainer,
	  "cbaumann_lowmass_CF.root");
	  */
  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("cbaumann_LMEEpp2010_CF",
		TList::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEpp2010_out.root");

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("cbaumann_EventStatPbPb2011",
		TH1D::Class(),
		AliAnalysisManager::kOutputContainer,
		"cbaumann_LMEEpp2010_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
