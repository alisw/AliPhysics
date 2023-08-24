AliAnalysisTask *AddTask_scheid_lowmass_r5(Bool_t getFromAlien=kFALSE,
					  TString cFileName = "Config_scheid_lowmass.C",
					  Char_t* outputFileName="LMEE.root",
					  ULong64_t triggerMask = AliVEvent::kINT7,
					  Bool_t pileupon = kFALSE,
					  Int_t wagon = 0
                                             )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_rbailhac_test", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hscheid/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if(hasMC) kMix = 0;

  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_rbailhac_lowmass")) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  //create task and add it to the manager (MB)
  TString appendix;
  if(wagon!=0) appendix += TString::Format("wagon%d_pileup%d",wagon,(Int_t)pileupon);
  else appendix += TString::Format("pileup%d",(Int_t)pileupon);
  printf("appendix %s\n", appendix.Data());
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_%s",appendix.Data()));
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);
  if(pileupon) task->SetRejectPileup();
  task->SetRandomizeDaughters(randomizeDau); //default kFALSE

  //Add event filter
  task->SetEventFilter( GetEventCuts(wagon) );

  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_rbailhac_lowmass(i);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);

  }//loop

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("tree_lowmass_%s", appendix.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("Histos_diel_lowmass_%s", appendix.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("CF_diel_lowmass_%s", appendix.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("scheid_lowmass_EventStat_%s", appendix.Data()),
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
