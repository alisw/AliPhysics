AliAnalysisTask *AddTask_lowmass_lowmasspp(){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lowmass", "No analysis manager found.");
    return 0;
  }


  Bool_t RunEMCtrigger = 0;
  Bool_t RunHighMulttrigger = 0;
  Bool_t RunMBtrigger = 1;

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);



  //Get the current train configuration
  
  //Directories for GSI train: 
  TString configBasePath("$TRAIN_ROOT/lowmass_LOWMASS/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  //Base Directory for GRID / LEGO Train
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";

  TString configFile("Config_lowmasspp.C");
  TString configFilePath(configBasePath+configFile);


  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
	  gROOT->LoadMacro(configFilePath.Data());


  if(RunMBtrigger){
  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *taskMB = new AliAnalysisTaskMultiDielectron("MultiDieMB");
  if (!hasMC) taskMB->UsePhysicsSelection();
  taskMB->SelectCollisionCandidates(AliVEvent::kMB);
  taskMB->SetRejectPileup();
	}

  if(RunHighMulttrigger){
  //create task and add it to the manager (HighMult)
  AliAnalysisTaskMultiDielectron *taskHighMult = new AliAnalysisTaskMultiDielectron("MultiDieHighMult");
  if (!hasMC) taskHighMult->UsePhysicsSelection();
  taskHighMult->SelectCollisionCandidates(AliVEvent::kHighMult);
  taskHighMult->SetRejectPileup();
	}

  if(RunEMCtrigger){
  //create task and add it to the manager (EMC1)
  AliAnalysisTaskMultiDielectron *taskEMC1 = new AliAnalysisTaskMultiDielectron("MultiDieEMC1");
  if (!hasMC) taskEMC1->UsePhysicsSelection();
  taskEMC1->SelectCollisionCandidates(AliVEvent::kEMC1);
  taskEMC1->SetRejectPileup();
}

          //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);

   if(RunMBtrigger)taskMB->SetEventFilter(eventCuts);
   if(RunHighMulttrigger) taskHighMult->SetEventFilter(eventCuts);
   if(RunEMCtrigger) taskEMC1->SetEventFilter(eventCuts);

    if(RunMBtrigger) mgr->AddTask(taskMB);
    if(RunHighMulttrigger)mgr->AddTask(taskHighMult);
    if(RunEMCtrigger)mgr->AddTask(taskEMC1);


  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
  if(RunMBtrigger){
    //MB
    AliDielectron *diel_lowMB = Config_lowmasspp(i);
    if(!diel_lowMB)continue;
    taskMB->AddDielectron(diel_lowMB);
	}
  if(RunHighMulttrigger){
    //HighMult
    AliDielectron *diel_lowHighMult = Config_lowmasspp(i);
    if(!diel_lowHighMult)continue;
    taskHighMult->AddDielectron(diel_lowHighMult);
}
	  if(RunEMCtrigger){
    //EMC1 
    AliDielectron *diel_lowEMC1 = Config_lowmasspp(i);
    if(!diel_lowEMC1)continue;
    taskEMC1->AddDielectron(diel_lowEMC1);
	}

  }//loop
  if(RunMBtrigger){
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmass.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmass.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("lowmass_lowmass_EventStat",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmass.root");

  mgr->ConnectInput(taskMB,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMB, 0, coutput1 );
  mgr->ConnectOutput(taskMB, 1, cOutputHist1);
  mgr->ConnectOutput(taskMB, 2, cOutputHist2);
  mgr->ConnectOutput(taskMB, 3, cOutputHist3);
  }
if(RunHighMulttrigger){
  //create output container (HighMult)
  AliAnalysisDataContainer *coutputHighMult1 =
    mgr->CreateContainer("tree_lowmassHighMult",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "default");
  
  AliAnalysisDataContainer *cOutputHistHighMult1 =
    mgr->CreateContainer("Histos_diel_lowmassHighMult",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassHighMult.root");

  AliAnalysisDataContainer *cOutputHistHighMult2 =
    mgr->CreateContainer("CF_diel_lowmassHighMult",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassHighMult.root");

  AliAnalysisDataContainer *cOutputHistHighMult3 =
    mgr->CreateContainer("lowmass_lowmass_EventStatHighMult",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassHighMult.root");
	}

	  if(RunEMCtrigger){

  //create output container (EMC1)
  AliAnalysisDataContainer *coutputEMC11 =
    mgr->CreateContainer("tree_lowmassEMC1",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "default");
  
  AliAnalysisDataContainer *cOutputHistEMC11 =
    mgr->CreateContainer("Histos_diel_lowmassEMC1",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassEMC1.root");

  AliAnalysisDataContainer *cOutputHistEMC12 =
    mgr->CreateContainer("CF_diel_lowmassEMC1",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassEMC1.root");

  AliAnalysisDataContainer *cOutputHistEMC13 =
    mgr->CreateContainer("lowmass_lowmass_EventStatEMC1",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "lowmass_lowmassEMC1.root");
}

//questionable regarding variable visibility
if(RunHighMulttrigger){
  mgr->ConnectInput(taskHighMult,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHighMult, 0, coutputHighMult1 );
  mgr->ConnectOutput(taskHighMult, 1, cOutputHistHighMult1);
  mgr->ConnectOutput(taskHighMult, 2, cOutputHistHighMult2);
  mgr->ConnectOutput(taskHighMult, 3, cOutputHistHighMult3);
}
  if(RunEMCtrigger){
  mgr->ConnectInput(taskEMC1,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEMC1, 0, coutputEMC11 );
  mgr->ConnectOutput(taskEMC1, 1, cOutputHistEMC11);
  mgr->ConnectOutput(taskEMC1, 2, cOutputHistEMC12);
  mgr->ConnectOutput(taskEMC1, 3, cOutputHistEMC13);
}


  return taskMB;

}
