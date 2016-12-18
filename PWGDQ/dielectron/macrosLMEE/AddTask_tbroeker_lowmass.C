AliAnalysisTask *AddTask_tbroeker_lowmass(Bool_t getFromAlien=kFALSE,
                                          TString cFileName = "Config_tbroeker_lowmass.C",
                                          Char_t* outputFileName="LMEE.root",
                                          ULong64_t triggerMask = AliVEvent::kINT7,
                                          Bool_t rejectPileup = kFALSE
                                         )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_tbroeker_lowmass", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/t/tbroker/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;
    
  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  gROOT->LoadMacro(configFilePath.Data());

  //Do we have an MC handler?
  hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if(hasMC) kMix = 0;
  
  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron("MultiDielectron");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);
  if(rejectPileup) task->SetRejectPileup(); //default kFALSE
  task->SetRandomizeDaughters(randomizeDau); //default kFALSE

  //Add event filter
  task->SetEventFilter( GetEventCuts() );
  
  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_tbroeker_lowmass(i);
    if(!diel_low)continue;
    task->AddDielectron(diel_low);

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
