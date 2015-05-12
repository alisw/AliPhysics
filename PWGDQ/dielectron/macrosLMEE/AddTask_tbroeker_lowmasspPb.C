AliAnalysisTask *AddTask_tbroeker_lowmasspPb(Bool_t getFromAlien=kFALSE,TString cFileName = "Config_lowmasspPb.C",Char_t* outputFileName="LMEEoutput.root"){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lowmass", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if(getFromAlien && (!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/t/tbroker/PWGDQ/dielectron/macrosLMEE/Config_lowmasspPb.C .")) ){
    TString configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
	  gROOT->LoadMacro(configFilePath.Data());

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *taskMB = new AliAnalysisTaskMultiDielectron("MultiDieMB");
  if (!hasMC) taskMB->UsePhysicsSelection();
//  taskMB->SelectCollisionCandidates(AliVEvent::kMB);
//  taskMB->SelectCollisionCandidates(AliVEvent::kINT7); //kINT7
  taskMB->SetTriggerMask(AliVEvent::kINT7);
//taskMB->SetRejectPileup();

  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);

  taskMB->SetEventFilter(eventCuts);
  mgr->AddTask(taskMB);


  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_lowMB = Config_lowmasspPb(i);
    if(!diel_lowMB)continue;
    diel_lowMB->SetNoPairing(kFALSE);
    taskMB->AddDielectron(diel_lowMB);

  }//loop

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("tbroeker_lowmass_EventStat",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  mgr->ConnectInput(taskMB,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMB, 0, coutput1 );
  mgr->ConnectOutput(taskMB, 1, cOutputHist1);
  mgr->ConnectOutput(taskMB, 2, cOutputHist2);
  mgr->ConnectOutput(taskMB, 3, cOutputHist3);
  


    return taskMB;

}
