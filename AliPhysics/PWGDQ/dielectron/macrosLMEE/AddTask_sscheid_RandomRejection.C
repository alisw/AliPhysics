AliAnalysisTask *AddTask_sscheid_RandomRejection(Bool_t getFromAlien=kFALSE,
                                                  Bool_t configsPreloaded=kFALSE,
                                                  TString cFileName = "Config_sscheid_lowmass.C",
                                                  Char_t* outputFileName="LMEE.root",
                                                  ULong64_t triggerMask = AliVEvent::kINT7
                                                 )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sscheid_RandomRejection", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if(getFromAlien && (!configsPreloaded) &&(!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hscheid/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  if (!configsPreloaded)
    gROOT->LoadMacro(configFilePath.Data());

  //create task and add it to the manager (MB)
  AliAnalysisTaskRandomRejection *task=new AliAnalysisTaskRandomRejection("RandomRejection");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);
//  taskMB->SetRejectPileup();
  task->SetRandomizeDaughters(randomizeDau); //default kFALSE

  //Add event filter
  task->SetEventFilter( GetEventCuts() );

  //task->SetPtFunc(PtFunc); // not working
  //printf(" task->GetPtFuncName(): %s\n",task->GetPtFuncName());
  task->SetPtExpr(RndmPtExpr);
  task->SetPtRange(RndmPtMin, RndmPtMax);
  task->SetEtaMax(RndmEtaMax);
  task->SetNTestpartPerEle(nTestpartPerEle);

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    //MB
    AliDielectron *diel_low = Config_sscheid_lowmass(i, kTRUE); //kTRUE -> "isRandomRejTask"
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("sscheid_RandomRejection_tree",
                       TTree::Class(),
                       AliAnalysisManager::kExchangeContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer("sscheid_RandomRejection_out",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer("sscheid_RandomRejection_CF",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer("sscheid_RandomRejection_EventStat",
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
