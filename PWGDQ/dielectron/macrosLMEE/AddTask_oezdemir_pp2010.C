AliAnalysisTask *AddTask_oezdemir_pp2010(Bool_t getFromAlien=kFALSE,TString cFileName = "Configpp2010Oezdemir.C",Char_t* outputFileName="LMEE.root"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_oezdemir_pp2010", "No analysis manager found.");
    return 0;
  }

//Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString configBasePath("$TRAIN_ROOT/oezdemir_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if (getFromAlien &&
      (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/mozdemir/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data())))
     ) {
        configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;
  
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //load dielectron configuration file
  gROOT->LoadMacro(configFilePath.Data());
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData");
  if (!hasMC) task->UsePhysicsSelection();

  //Add event filter
  task->SetEventFilter( GetEventCuts() );

  //Min Bias?
  task->SetTriggerMask(AliVEvent::kMB);
  task->SetRandomizeDaughters(kTRUE);
  mgr->AddTask(task);

  //If MC available decide which pdg codes are tested:
  for (Int_t i=0; i<nDie; ++i){
  AliDielectron *lowmass=Configpp2010(i,kFALSE); //kFALSE->RandomRjection off
  task->AddDielectron(lowmass);
  printf("add: %s\n",lowmass->GetName());
  }
	//create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("oezdemir_pp2010_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("oezdemir_pp2010_out",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("oezdemir_pp2010_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("oezdemir_pp2010_EventStat",
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
