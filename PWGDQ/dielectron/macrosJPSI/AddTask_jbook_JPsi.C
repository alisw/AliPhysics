AliAnalysisTask *AddTask_jbook_JPsi(TString config="1",
						TString cfg="ConfigJpsi_jb_PbPb.C",
				    Bool_t gridconf=kFALSE,
				    Bool_t hasMC=kFALSE,
				    ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB,
						Bool_t bMultiToSingle=kTRUE){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jbook_JPsi", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  TString list = gSystem->Getenv("LIST");
  if(!list.IsNull()) {
    if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
    if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;
  }

  //Do we have an AOD handler?
  Bool_t isAOD=(mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class() ? kTRUE : kFALSE);

  // set AOD debug levels
  if(isAOD) {
    mgr->AddClassDebug("AliAODTrack", AliLog::kFatal);
    mgr->AddClassDebug("AliAODpidUtil", AliLog::kInfo);
  }

  //set config file name
  TString configFile("");
  printf("%s \n",gSystem->pwd());
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsi_jb_PbPb.C";

  // the different paths
  TString gsiPath("$TRAIN_ROOT/jbook_jpsi/");
  TString alienPath("alien:///alice/cern.ch/user/j/jbook/PWGDQ/dielectron/macrosJPSI/");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/");

  ////////// >>>>>>>>>> gsi config
  if (!trainRoot.IsNull())  configFile=gsiPath.Data();
  ////////// >>>>>>>>>> alien config
  else if(!gSystem->Exec(Form("alien_cp %s/%s .",alienPath.Data(),cfg.Data()))) {
    gSystem->Exec(Form("ls -l %s",gSystem->pwd()));
    configFile=gSystem->pwd();
  }
  else {
    printf("ERROR: couldn't copy file %s/%s from grid \n", alienPath.Data(),cfg.Data() );
    return;
  }
  ///////// >>>>>>>>> aliroot config
  if(!gridconf && trainRoot.IsNull()) configFile=alirootPath.Data();
  ///////// add config to path
  configFile+="/";
  configFile+=cfg.Data();

  // trigger selection
  ULong64_t triggerSets[]={AliVEvent::kCentral , AliVEvent::kSemiCentral , AliVEvent::kMB,
			   AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB};
  const char* triggerNames[]={"Central","SemiCentral","MB","ALL"};
  const char* onlineRejection[]={"","CCENT","",""};

  // find out the configured triggers
  Int_t j=0;
  for(j=0; j<4; j++) {
    if(triggers!=triggerSets[j]) continue;
    else break;
  }

  // print overall configuration
  printf("production: %s MC: %d \n",  list.Data(),hasMC);
  printf("triggers:   %s \n",         triggerNames[j]  );
  printf("config:     %s Grid: %d \n",configFile.Data(),gridconf);

  //create task(s)
  AliAnalysisTaskMultiDielectron *task;
  if(!bMultiToSingle) {
    // create one multi task
    task = new AliAnalysisTaskMultiDielectron(Form("MultiDieJB"));
    task->SetBeamEnergy(1380.);
    task->SetTriggerMask(triggers);
    if(strlen(onlineRejection[j])) task->SetFiredTriggerName(onlineRejection[j],kTRUE);
    if(!hasMC) task->UsePhysicsSelection();
  }

  // event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("vertex","vertex");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  if(hasMC) eventCuts->SetVertexZ(-10.,+10.); //for data this is done by in the config
  eventCuts->SetCentralityRange(0,90.);
  eventCuts->Print();
  if(!bMultiToSingle)       task->SetEventFilter(eventCuts);

  //load dielectron configuration file (only once)
  TString checkconfig="ConfigJpsi_jb_PbPb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  //define default output container
  TString containerName = "JPSI.root";

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i) { //nDie defined in config file

    //only configs switched ON will pass
    if(config.Length()<=i || config(i,1)!="1") { printf(" %d switched OFF \n",i); continue; }

    // load configuration
    AliDielectron *jpsi=ConfigJpsi_jb_PbPb(i,hasMC,triggers);
    if(!jpsi) continue;

    // create unique title
    TString unitit = Form("%s_%s",triggerNames[j],jpsi->GetName());

    // create single tasks instead of one multi task (decreasing size of CF container)
    if(bMultiToSingle) {
      task = new AliAnalysisTaskMultiDielectron(Form("MultiDieJB_%s",unitit.Data()));
      task->SetBeamEnergy(1380.);
      task->SetTriggerMask(triggers);
      if(strlen(onlineRejection[j])) task->SetFiredTriggerName(onlineRejection[j],kTRUE);
      if(!hasMC) task->UsePhysicsSelection();
    }

    // add dielectron to the task and manager
    task->AddDielectron(jpsi);

    // multiple output connection
    if(bMultiToSingle) {
      task->SetEventFilter(eventCuts);
      mgr->AddTask(task);

      //create output sub containers
      unitit.Prepend("jbook_QA_");
      AliAnalysisDataContainer *cOutputHist1 =
	mgr->CreateContainer(unitit.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());
      unitit.ReplaceAll("_QA_","_CF_");
      AliAnalysisDataContainer *cOutputHist2 =
	mgr->CreateContainer(unitit.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());
      unitit.ReplaceAll("_CF_","_EventStat_");
      AliAnalysisDataContainer *cOutputHist3 =
	mgr->CreateContainer(unitit.Data(), TH1D::Class(), AliAnalysisManager::kOutputContainer,containerName.Data());

      mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      //  mgr->ConnectOutput(task, 0, coutput1 );
      mgr->ConnectOutput(task, 1, cOutputHist1);
      mgr->ConnectOutput(task, 2, cOutputHist2);
      mgr->ConnectOutput(task, 3, cOutputHist3);
    }

    printf(" %s added\n",jpsi->GetName());

  } //end : loop over configs

  
  // multiple output connection
  if(!bMultiToSingle) {
    mgr->AddTask(task);
    
    //create output sub containers
    AliAnalysisDataContainer *cOutputHist1 =
      mgr->CreateContainer("jbook_QA", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());
    AliAnalysisDataContainer *cOutputHist2 =
      mgr->CreateContainer("jbook_CF", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());
    AliAnalysisDataContainer *cOutputHist3 =
      mgr->CreateContainer("jbook_EventStat", TH1D::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());

    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      //  mgr->ConnectOutput(task, 0, coutput1 );
    mgr->ConnectOutput(task, 1, cOutputHist1);
    mgr->ConnectOutput(task, 2, cOutputHist2);
    mgr->ConnectOutput(task, 3, cOutputHist3);
  }


  return task;
}
