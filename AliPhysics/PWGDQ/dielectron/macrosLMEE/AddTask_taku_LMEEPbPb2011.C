/*
  in real data. argument should be kFALSE, kFALSE, kFALSE,
  in siml data. argument should be kFALSE, kTRUE, kTRUE,

*/
AliAnalysisTask *AddTask_taku_LMEEPbPb2011(Bool_t runRejection=kFALSE, Bool_t setMC=kFALSE,Bool_t enableCF=kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_taku_LMEEPbPb2011", "No analysis manager found.");
    return 0;
  }

  //create config File names: TRAIN_ROOT is for running on GSI train, 
  // ALICE_ROOT for CERN Lego trains
  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");                                                                            
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";
  TString configFile("ConfigTakuLMEEPbPb2011.C");
  TString configLMEECutLib("LMEECutLibTaku.C");

  TString configFilePath(configBasePath+configFile);
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);

  //AOD Usage currently tested with separate task, to be merged
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }

  //Do we have an MC handler?
  Bool_t hasMC=setMC;
  if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0){
    hasMC=kTRUE;
  }
  
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());

  LMEECutLib* cutlib = new LMEECutLib();
  cutlib->SetMCFlag(hasMC);
  AliAnalysisTaskMultiDielectronTG *task=new AliAnalysisTaskMultiDielectronTG("MultiDiEData");

  ////default cutter defined in ConfigTakuLMEEPbPb2011.C
  Int_t PairCutTypeDef[20]={
			    0,0,0, //no pair cuts
			    1,1,1, //reject from arrays by op cuts
			    2,2,2, //reject from arrays by phiv cuts
			    3,3,3, //pair-by-pair cuts by op
			    4,4,4, //pair-by-pair cuts by phiv
			    0,0,0,
			    0,0};

  if (!hasMC){ 
    task->UsePhysicsSelection();
  }
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kPbPb2011TPCandTOF)); //
  task->RejectConversion(2.0);
  task->RejectOP(0.035);
  task->SetPairCuts(PairCutTypeDef);
  task->EnableV0mixing(kFALSE);

  //load dielectron configuration file
  //add dielectron analysis with different cuts to the task
  if (runRejection) {
    AliDielectron *lowmass4=ConfigTakuLMEEPbPb2011(4,hasMC,enableCF);
    task->AddDielectron(lowmass4);
    printf("add: %s\n",lowmass4->GetName());
  }
  else {
    //////// this is for test 
    AliDielectron *lowmass1=ConfigTakuLMEEPbPb2011(1,hasMC,enableCF);
    lowmass1->SetUseKF(kFALSE);
    task->AddDielectron(lowmass1);
    printf("add: %s\n",lowmass1->GetName());

    AliDielectron *lowmass2=ConfigTakuLMEEPbPb2011(2,hasMC,enableCF);
    lowmass2->SetUseKF(kFALSE);
    task->AddDielectron(lowmass2);
    printf("add: %s\n",lowmass2->GetName());

    AliDielectron *lowmass3=ConfigTakuLMEEPbPb2011(3,hasMC,enableCF);
    lowmass3->SetUseKF(kFALSE);
    task->AddDielectron(lowmass3);
    printf("add: %s\n",lowmass3->GetName());

    ///////////////////////////

    AliDielectron *lowmass4=ConfigTakuLMEEPbPb2011(4,hasMC,enableCF);
    lowmass4->SetUseKF(kFALSE);
    task->AddDielectron(lowmass4);
    printf("add: %s\n",lowmass4->GetName());

    AliDielectron *lowmass5=ConfigTakuLMEEPbPb2011(5,hasMC,enableCF);
    lowmass5->SetUseKF(kFALSE);
    task->AddDielectron(lowmass5);
    printf("add: %s\n",lowmass5->GetName());

    AliDielectron *lowmass6=ConfigTakuLMEEPbPb2011(6,hasMC,enableCF);
    lowmass6->SetUseKF(kFALSE);
    task->AddDielectron(lowmass6);
    printf("add: %s\n",lowmass6->GetName());


    ///////////////////////////

    AliDielectron *lowmass7=ConfigTakuLMEEPbPb2011(7,hasMC,enableCF);
    lowmass7->SetUseKF(kFALSE);
    task->AddDielectron(lowmass7);
    printf("add: %s\n",lowmass7->GetName());

    AliDielectron *lowmass8=ConfigTakuLMEEPbPb2011(8,hasMC,enableCF);
    lowmass8->SetUseKF(kFALSE);
    task->AddDielectron(lowmass8);
    printf("add: %s\n",lowmass8->GetName());

    AliDielectron *lowmass9=ConfigTakuLMEEPbPb2011(9,hasMC,enableCF);
    lowmass9->SetUseKF(kFALSE);
    task->AddDielectron(lowmass9);
    printf("add: %s\n",lowmass9->GetName());

    ///////////////////////////

    AliDielectron *lowmass10=ConfigTakuLMEEPbPb2011(10,hasMC,enableCF);
    lowmass10->SetUseKF(kFALSE);
    task->AddDielectron(lowmass10);
    printf("add: %s\n",lowmass10->GetName());

    AliDielectron *lowmass11=ConfigTakuLMEEPbPb2011(11,hasMC,enableCF);
    lowmass11->SetUseKF(kFALSE);
    task->AddDielectron(lowmass11);
    printf("add: %s\n",lowmass11->GetName());

    AliDielectron *lowmass12=ConfigTakuLMEEPbPb2011(12,hasMC,enableCF);
    lowmass12->SetUseKF(kFALSE);
    task->AddDielectron(lowmass12);
    printf("add: %s\n",lowmass12->GetName());

    ///////////////////////////


    AliDielectron *lowmass13=ConfigTakuLMEEPbPb2011(13,hasMC,enableCF);
    lowmass13->SetUseKF(kFALSE);
    task->AddDielectron(lowmass13);
    printf("add: %s\n",lowmass13->GetName());

    AliDielectron *lowmass14=ConfigTakuLMEEPbPb2011(14,hasMC,enableCF);
    lowmass14->SetUseKF(kFALSE);
    task->AddDielectron(lowmass14);
    printf("add: %s\n",lowmass14->GetName());

    AliDielectron *lowmass15=ConfigTakuLMEEPbPb2011(15,hasMC,enableCF);
    lowmass15->SetUseKF(kFALSE);
    task->AddDielectron(lowmass15);
    printf("add: %s\n",lowmass15->GetName());



  }

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("taku_LMEEPbPb2011_tree",
	TTree::Class(),
	AliAnalysisManager::kExchangeContainer,
	"taku_LMEEPbPb2011_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("taku_LMEEPbPb2011_out",
	TList::Class(),
	AliAnalysisManager::kOutputContainer,
	"taku_LMEEPbPb2011_out.root");
  
  AliAnalysisDataContainer *cOutputHist2 = 0x0;
  if (enableCF) {
    cOutputHist2 = 
      mgr->CreateContainer("taku_LMEEPbPb2011_CF",
	  TList::Class(),
	  AliAnalysisManager::kOutputContainer,
	  "taku_LMEEPbPb2011_out.root");

  }
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("taku_EventStatPbPb2011",
	TH1D::Class(),
	AliAnalysisManager::kOutputContainer,
	"taku_LMEEPbPb2011_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  if (enableCF) {
    mgr->ConnectOutput(task, 2, cOutputHist2);
  }
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
