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
  TString configLMEECutLib("LMEECutLib.C");

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
    enableCF=kTRUE;
  }
  


  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data()))
    gROOT->LoadMacro(configLMEECutLibPath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());

  LMEECutLib* cutlib = new LMEECutLib();
  cutlib->SetMCFlag(hasMC);
  AliAnalysisTaskMultiDielectronTG *task=new AliAnalysisTaskMultiDielectronTG("MultiDiEData");

  if (!hasMC){ 
    task->UsePhysicsSelection();
  }
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  task->SetEventFilter(cutlib->GetEventCuts(LMEECutLib::kPbPb2011TPCandTOF)); //
  task->reject_conversion(2.0);
  task->enable_v0mixing(kFALSE);

  //load dielectron configuration file
  //add dielectron analysis with different cuts to the task
  if (runRejection) {
    AliDielectron *lowmass4=ConfigTakuLMEEPbPb2011(4,hasMC,enableCF);
    task->AddDielectron(lowmass4);
    printf("add: %s\n",lowmass4->GetName());
  }
  else {
    AliDielectron *lowmass3=ConfigTakuLMEEPbPb2011(3,hasMC,enableCF);
    task->AddDielectron(lowmass3);
    printf("add: %s\n",lowmass3->GetName());

    AliDielectron *lowmass1=ConfigTakuLMEEPbPb2011(1,hasMC,enableCF);
    task->AddDielectron(lowmass1);
    printf("add: %s\n",lowmass1->GetName());

    AliDielectron *lowmass7=ConfigTakuLMEEPbPb2011(7,hasMC,enableCF);
    task->AddDielectron(lowmass7);
    printf("add: %s\n",lowmass7->GetName());

    AliDielectron *lowmass2=ConfigTakuLMEEPbPb2011(2,hasMC,enableCF);
    task->AddDielectron(lowmass2);
    printf("add: %s\n",lowmass2->GetName());
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
