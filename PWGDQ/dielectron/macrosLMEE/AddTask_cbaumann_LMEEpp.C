AliAnalysisTask *AddTask_cbaumann_LMEEpp(Bool_t enablePS=kTRUE){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_cbaumann", "No analysis manager found.");
    return 0;
  }


  Bool_t RunEMCtrigger = 0;
  Bool_t RunHighMulttrigger = 0;
  Bool_t RunMBtrigger = 1;

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);



  //Get the current train configuration
  
  //Directories for GSI train: 
  TString configBasePath("$TRAIN_ROOT/cbaumann_dielectron/");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");

  //Base Directory for GRID / LEGO Train
  if (trainRoot.IsNull()) configBasePath= "$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/";

  TString configFile("Config_lowmasspp.C");
  TString configFilePath(configBasePath+configFile);


  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
	  gROOT->LoadMacro(configFilePath.Data());


  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *taskMB = new AliAnalysisTaskMultiDielectron("MultiDieMB");
  if ((!hasMC) && enablePS) taskMB->UsePhysicsSelection();
  //taskMB->SelectCollisionCandidates(AliVEvent::kMB);
  taskMB->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kMB);
  taskMB->SelectCollisionCandidates(AliVEvent::kAny);

//  taskMB->SetRejectPileup();

          //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);

  taskMB->SetEventFilter(eventCuts);

   mgr->AddTask(taskMB);


  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_lowMB = Config_lowmasspp(i);
    if(!diel_lowMB)continue;
    taskMB->AddDielectron(diel_lowMB);

  }//loop
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_cb_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_cbaumann_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_cbaumann_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("cbaumann_lowmass_EventStat",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  mgr->ConnectInput(taskMB,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMB, 0, coutput1 );
  mgr->ConnectOutput(taskMB, 1, cOutputHist1);
  mgr->ConnectOutput(taskMB, 2, cOutputHist2);
  mgr->ConnectOutput(taskMB, 3, cOutputHist3);

  return taskMB;

}
