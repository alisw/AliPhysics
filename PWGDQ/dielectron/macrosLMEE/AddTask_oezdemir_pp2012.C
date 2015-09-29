AliAnalysisTask *AddTask_oezdemir_pp2012(Bool_t getFromAlien=kFALSE,Char_t* outputFileName="LMEE.root"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_oezdemir_pp2012", "No analysis manager found.");
    return 0;
  }

//Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString configBasePath("$TRAIN_ROOT/oezdemir_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien &&
      (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/m/mozdemir/PWGDQ/dielectron/macrosLMEE/Configpp2012Oezdemir.C"))
     ) {
        configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFile("Configpp2012Oezdemir.C");

  TString configFilePath(configBasePath+configFile);

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
  //2012 Min Bias?
  task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kMB+AliVEvent::kINT8);

  mgr->AddTask(task);


  //load dielectron configuration file
  gROOT->LoadMacro(configFilePath.Data());

  //If MC available decide which pdg codes are tested:

  AliDielectron *lowmass0=Configpp2012Oezdemir(0,hasMC);
  task->AddDielectron(lowmass0);
  printf("add: %s\n",lowmass0->GetName());

	//create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("oezdemir_pp2012_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("oezdemir_pp2012_out",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("oezdemir_pp2012_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("oezdemir_pp2012_EventStat",
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
