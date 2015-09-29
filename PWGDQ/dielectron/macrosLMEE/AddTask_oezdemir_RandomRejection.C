AliAnalysisTask *AddTask_oezdemir_RandomRejection(Bool_t getFromAlien=kFALSE, Bool_t configsPreloaded=kFALSE, TString cFileName = "Configpp2010Oezdemir.C",Char_t* outputFileName="LMEE.root"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_oezdemir_RandomRejection", "No analysis manager found.");
    return 0;
  }

//Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString configBasePath("$TRAIN_ROOT/oezdemir_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien && !configsPreloaded &&
      (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/mozdemir/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data())))
      ) {
    configBasePath=Form("%s/",gSystem->pwd());
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
 
  //load dielectron configuration file
  if (!configsPreloaded) {
    TString configFilePath(configBasePath+cFileName);
    std::cout << "loading config file:  " << configFilePath << std::endl;
    gROOT->LoadMacro(configFilePath.Data());
  } else {
    std::cout << "using preloaded config file." << std::endl;
  }
   
  //create task and add it to the manager
  AliAnalysisTaskRandomRejection *task=new AliAnalysisTaskRandomRejection("MultiDiE_RandomRejection");
  if (!hasMC) task->UsePhysicsSelection();

  //Add event filter
  task->SetEventFilter( GetEventCuts() );  
  //AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  //eventCuts->SetRequireVertex();
  //eventCuts->SetVertexZ(-10.,10.);
  //eventCuts->SetMinVtxContributors(1); 
  
  //task->SetEventFilter(eventCuts);

  //Min Bias?
  task->SetTriggerMask(AliVEvent::kMB);
  task->SetRandomizeDaughters(kFALSE);
  mgr->AddTask(task);


  for (Int_t i=0; i<nDie; ++i){
  AliDielectron *lowmass=Configpp2010(i,kTRUE); //kTRUE->Random Rejection on
  task->AddDielectron(lowmass);
  printf("add: %s\n",lowmass->GetName());
  }
	//create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("oezdemir_RandomRejection_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("oezdemir_RandomRejection_out",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("oezdemir_RandomRejection_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("oezdemir_RandomRejection_EventStat",
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
