AliAnalysisTask *AddTask_jjung_lowmass(Bool_t getFromAlien=kFALSE,
                                      TString cFileName = "Config_jjung_lowmass.C",
                                      Char_t* outputFileName="LMEE.root",
                                      Int_t wagonnr=0                                    
                                      )
{
  std::cout << "AddTask_jjung_lowmass" << std::endl;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jjung_lowmass", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;
  
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  gROOT->LoadMacro(configFilePath.Data());

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_%d", wagonnr));
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);

  if(nSetV0ANDTrigger)
    task->SetTriggerOnV0AND(kTRUE); // only for cross-check
  
  // SPD pile-up rejection with 3 contributors is applied only for min.bias analysis
  if(nSetPileupRejection)
    task->SetRejectPileup(kTRUE); // rejectPileup == kTRUE by default for min.bias

  task->SetRandomizeDaughters(randomizeDau); //default kFALSE

  //Add event filter
  task->SetEventFilter( GetEventCuts() );  
  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_jjung_lowmass(i,kFALSE,wagonnr);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
  }//loop

  //create output container
  AliAnalysisDataContainer *coutput1 = 0x0;
  AliAnalysisDataContainer *cOutputHist1 = 0x0;
  AliAnalysisDataContainer *cOutputHist2 = 0x0;
  AliAnalysisDataContainer *cOutputHist3 = 0x0;


  if(wagonnr == 0){
    coutput1 = mgr->CreateContainer("tree_lowmass", TTree::Class(),
                         AliAnalysisManager::kExchangeContainer, outputFileName);
  
    cOutputHist1 = mgr->CreateContainer("Histos_diel_lowmass", TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);

    cOutputHist2 = mgr->CreateContainer("CF_diel_lowmass", TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);

    cOutputHist3 = mgr->CreateContainer("jjung_lowmass_EventStat", TH1D::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);
  }

  else{
    coutput1 = mgr->CreateContainer(Form("tree_lowmass_%d", wagonnr), TTree::Class(),
                         AliAnalysisManager::kExchangeContainer, outputFileName);
  
    cOutputHist1 = mgr->CreateContainer(Form("Histos_diel_lowmass_%d", wagonnr), TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);

    cOutputHist2 = mgr->CreateContainer(Form("CF_diel_lowmass_%d", wagonnr), TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);

    cOutputHist3 = mgr->CreateContainer(Form("jjung_lowmass_EventStat_%d", wagonnr), TH1D::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName);

  }

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
  
}

