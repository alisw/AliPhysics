AliAnalysisTask *AddTask_oezdemir_LOWMASS(Bool_t getFromAlien=kFALSE){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_oezdemir_LOWMASS", "No analysis manager found.");
    return 0;
  }

//Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString configBasePath("$TRAIN_ROOT/oezdemir_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cbaumann/PWGDQ/dielectron/macrosLMEE/ConfigLowMassDiEOezdemir.C"))
     ) {
        configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFile("ConfigLowMassDiEOezdemir.C");

  TString configFilePath(configBasePath+configFile);
/*
  //AOD Usage currently tested with separate task, to be merged
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
        ::Info("AddTaskLMEEPbPb2011", "no dedicated AOD configuration");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
        ::Info("AddTaskLMEEPbPb2011AOD","switching on ESD specific code");
        bESDANA=kTRUE;
  }
*/
/*
//AOD Usage currently not allows-----------------------------------------

  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskJPSI", "Using AOD configuration");
  //  configFile="$TRAIN_ROOT/util/dielectron/dielectron/macros/ConfigJpsi2eeDataAOC.C";
  }
*/
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();

//Add event filter
AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
eventCuts->SetRequireVertex();
eventCuts->SetVertexZ(-10.,10.);
eventCuts->SetMinVtxContributors(1);

  task->SetEventFilter(eventCuts);
  //2010/2011	Min Bias?
  task->SetTriggerMask(AliVEvent::kMB);

  mgr->AddTask(task);


  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());

  //If MC available decide which pdg codes are tested:

  AliDielectron *lowmass0=ConfigLowMassDiEOezdemir(0,hasMC);
  task->AddDielectron(lowmass0);
  printf("add: %s\n",lowmass0->GetName());
/*
  AliDielectron *lowmass1=ConfigLowMassDiEOezdemir(1,hasMC);
  task->AddDielectron(lowmass1);
  printf("add: %s\n",lowmass1->GetName());


  AliDielectron *lowmass2=ConfigLowMassDiE(2,hasMC);
  task->AddDielectron(lowmass2);
  printf("add: %s\n",lowmass2->GetName());

  AliDielectron *lowmass3=ConfigLowMassDiE(3,hasMC);
  task->AddDielectron(lowmass3);
  printf("add: %s\n",lowmass3->GetName());

  AliDielectron *lowmass4=ConfigLowMassDiE(4,hasMC);
  task->AddDielectron(lowmass4);
  printf("add: %s\n",lowmass4->GetName());

  AliDielectron *lowmass5=ConfigLowMassDiE(5,hasMC);
  task->AddDielectron(lowmass5);
  printf("add: %s\n",lowmass5->GetName());

  AliDielectron *lowmass6=ConfigLowMassDiE(6,hasMC);
  task->AddDielectron(lowmass6);
  printf("add: %s\n",lowmass6->GetName());

  AliDielectron *lowmass7=ConfigLowMassDiE(7,hasMC);
  task->AddDielectron(lowmass7);
  printf("add: %s\n",lowmass7->GetName());
*/
	//create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("oezdemir_LOWMASS_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "oezdemir_LOWMASS_default.root");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("oezdemir_LOWMASS_out",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "oezdemir_LOWMASS_out.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("oezdemir_LOWMASS_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "oezdemir_LOWMASS_out.root");
//                         "oezdemir_LOWMASS_CF.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("oezdemir_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "oezdemir_LOWMASS_out.root");

  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
