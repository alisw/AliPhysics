AliAnalysisTask *AddTask_jjung_RandomRejection(Bool_t getFromAlien=kFALSE,
                                                  TString cFileName = "Config_jjung_lowmass.C",
                                                  Char_t* outputFileName="LMEE.root",
                                                 )
{ 

  std::cout << "AddTask_jjung_RandomRejection" << std::endl;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jjung_RandomRejection", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "/data4/jung/localLegotrainRandRejection/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);

  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
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
    AliDielectron *diel_low = Config_jjung_lowmass(i, kTRUE); //kTRUE -> "isRandomRejTask"
    if(!diel_low)continue;
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }//loop

  mgr->AddTask(task);

  //create output container
  AliAnalysisDataContainer *coutput1 =
        mgr->CreateContainer("jjung_RandomRejection_tree",
                       TTree::Class(),
                       AliAnalysisManager::kExchangeContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
        mgr->CreateContainer("jjung_RandomRejection_out",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist2 =
        mgr->CreateContainer("jjung_RandomRejection_CF",
                       TList::Class(),
                       AliAnalysisManager::kOutputContainer,
                       outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
        mgr->CreateContainer("jjung_RandomRejection_EventStat",
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
