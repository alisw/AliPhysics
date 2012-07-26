/*
  in real data. argument should be kFALSE, kFALSE, kFALSE,
  in siml data. argument should be kFALSE, kTRUE, kTRUE,

*/
AliAnalysisTask *AddTask_taku_LMEEPbPb2011SemiCent2(Bool_t runRejection=kFALSE, Bool_t setMC=kFALSE,Bool_t enableCF=kFALSE){
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
  
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())){
    if(!gROOT->GetClass("LMEECutLibTaku")){
      gROOT->LoadMacro(configLMEECutLibPath.Data());
    }
  }
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data())){
    gROOT->LoadMacro(configFilePath.Data());
  }
  LMEECutLibTaku* cutlib = new LMEECutLibTaku();
  cutlib->SetMCFlag(hasMC);
  AliAnalysisTaskMultiDielectronTG *task=new AliAnalysisTaskMultiDielectronTG("MultiDiETGDataSemiCent2");

  ////default cutter defined in ConfigTakuLMEEPbPb2011.C
  Int_t PairCutTypeDef[20]={0,
			    0,0,0, //no pair cuts
			    1,1,1, //reject from arrays by op cuts
			    2,2,2, //reject from arrays by phiv cuts
			    3,3,3, //pair-by-pair cuts by op
			    4,4,4, //pair-by-pair cuts by phiv
			    0,0,0,
			    0};

  Int_t PairCutType[20]={0};
  PairCutType[0] = PairCutTypeDef[3];
  PairCutType[1] = PairCutTypeDef[6];
  PairCutType[2] = PairCutTypeDef[9];
  PairCutType[3] = PairCutTypeDef[12];
  PairCutType[4] = PairCutTypeDef[15];

  if (!hasMC){ 
    task->UsePhysicsSelection();
  }
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLibTaku::kPbPb2011TPCandTOF)); //
  task->RejectConversion(2.0, 0.3);
  task->RejectOP(0.035);
  task->SetPairCuts(PairCutType);
  task->EnableV0mixing(kFALSE);
  task->SetRejBGPairs(kFALSE, kTRUE);

  //load dielectron configuration file
  //add dielectron analysis with different cuts to the task
  if (runRejection) {
    AliDielectron *lowmass4=ConfigTakuLMEEPbPb2011(4,hasMC,enableCF);
    task->AddDielectron(lowmass4);
    printf("add: %s\n",lowmass4->GetName());
  }
  else {
    //////// this is for test 
    AliDielectron *lowmass3=ConfigTakuLMEEPbPb2011(3,hasMC,enableCF);
    lowmass3->SetUseKF(kFALSE);
    task->AddDielectron(lowmass3);
    printf("add: %s\n",lowmass3->GetName());

    ///////////////////////////

    AliDielectron *lowmass6=ConfigTakuLMEEPbPb2011(6,hasMC,enableCF);
    lowmass6->SetUseKF(kFALSE);
    task->AddDielectron(lowmass6);
    printf("add: %s\n",lowmass6->GetName());


    ///////////////////////////


    AliDielectron *lowmass9=ConfigTakuLMEEPbPb2011(9,hasMC,enableCF);
    lowmass9->SetUseKF(kFALSE);
    task->AddDielectron(lowmass9);
    printf("add: %s\n",lowmass9->GetName());


    ///////////////////////////

    AliDielectron *lowmass12=ConfigTakuLMEEPbPb2011(12,hasMC,enableCF);
    lowmass12->SetUseKF(kFALSE);
    task->AddDielectron(lowmass12);
    printf("add: %s\n",lowmass12->GetName());

    ///////////////////////////



    AliDielectron *lowmass15=ConfigTakuLMEEPbPb2011(15,hasMC,enableCF);
    lowmass15->SetUseKF(kFALSE);
    task->AddDielectron(lowmass15);
    printf("add: %s\n",lowmass15->GetName());

  }

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("taku_LMEEPbPb2011_semicent2_tree",
	TTree::Class(),
	AliAnalysisManager::kExchangeContainer,
	"taku_LMEEPbPb2011_semicent2_default.root");

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("taku_LMEEPbPb2011_semicent2_out",
	TList::Class(),
	AliAnalysisManager::kOutputContainer,
	"taku_LMEEPbPb2011_semicent2_out.root");
  
  AliAnalysisDataContainer *cOutputHist2 = 0x0;
  if (enableCF) {
    cOutputHist2 = 
      mgr->CreateContainer("taku_LMEEPbPb2011_semicent2_CF",
	  TList::Class(),
	  AliAnalysisManager::kOutputContainer,
	  "taku_LMEEPbPb2011_semicent2_out.root");

  }
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("taku_EventStatPbPb2011_semicent2",
	TH1D::Class(),
	AliAnalysisManager::kOutputContainer,
	"taku_LMEEPbPb2011_semicent2_out.root");


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  if (enableCF) {
    mgr->ConnectOutput(task, 2, cOutputHist2);
  }
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
