AliAnalysisTask *AddTask_pdill_JPsi(TString config="1",
				    TString cfg="ConfigJpsi_pd_PbPb.C",
				    Bool_t gridconf=kFALSE,
				    Bool_t hasMC=kFALSE,
				    ULong64_t triggers=AliVEvent::kAnyINT,
						TString period=""
				    ){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_pdill_JPsi", "No analysis manager found.");
    return 0;
  }
  printf("------------------------------------------------\n");

  //Do we have an MC handler?
  TString list = gSystem->Getenv("REFERENCE_PRODUCTION");
  if(!list.IsNull()) {
    if( list.Contains("LHC16g1") || list.Contains("LHC16j1") || list.Contains("LHC16f5") ) hasMC=kTRUE;
  }
	TString runNumString = gSystem->Getenv("RUNNO");
	Int_t runnumber = runNumString.Atoi();

	printf("HalloWelt REF_PROD: %s Run: %d %d\n", list.Data(), runnumber, __LINE__);

  //Do we have an AOD handler?
  Bool_t isAOD=(mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class() ? kTRUE : kFALSE);

  // set AOD debug levels
  if(isAOD) {
    mgr->AddClassDebug("AliAODTrack", AliLog::kFatal);
    mgr->AddClassDebug("AliAODpidUtil", AliLog::kInfo);
  }

  //set config file name
  TString configFile("");
  printf("pwd:        %s \n",gSystem->pwd());
  if(cfg.IsNull()) cfg="ConfigJpsi_pd_pp.C";

  // the different paths
  TString alienPath("alien:///alice/cern.ch/user/p/pdillens/PWGDQ/dielectron/macrosJPSI/");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/");

  ////////// >>>>>>>>>> alien config
  if(gridconf) {
    if(!gSystem->Exec(Form("alien_cp %s/%s .",alienPath.Data(),cfg.Data()))) {
      gSystem->Exec(Form("ls -l %s",gSystem->pwd()));
      configFile=gSystem->pwd();
    }
    else {
      printf("ERROR: couldn't copy file %s/%s from grid \n", alienPath.Data(),cfg.Data() );
      return;
    }
  }
  ///////// >>>>>>>>> aliroot config
  else if(!gridconf) configFile=alirootPath.Data();
  ///////// add config to path
  configFile+="/";
  configFile+=cfg.Data();

  // trigger selection
  ULong64_t triggerSets[]={AliVEvent::kAnyINT};
  const char* triggerNames[]={"MinBias"};
  Int_t j=0;

  // print overall configuration
  printf("production: %s MC: %d \n",  list.Data(),hasMC);
  printf("triggers:   %s \n",         triggerNames[j]  );
  printf("config:     %s Grid: %d \n",configFile.Data(),gridconf);
  printf("------------------------------------------------\n");

  //create task(s)
  AliAnalysisTaskMultiDielectron *task;
  // create one multi task
  task = new AliAnalysisTaskMultiDielectron(Form("MultiDieJB"));
  task->SetBeamEnergy(2510.);
  task->SetTriggerMask(triggers);
  if(!hasMC) task->UsePhysicsSelection();

  // event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("vertex","vertex");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  if(hasMC) eventCuts->SetVertexZ(-10.,+10.); //for data this is done by in the config
  eventCuts->Print();
  task->SetEventFilter(eventCuts);

  //load dielectron configuration file (only once)
  TString checkconfig="ConfigJpsi_pd_pp";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  //define default output container
  TString containerName = "JPSI.root";

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i) { //nDie defined in config file

    //only configs switched ON will pass
    if(config.Length()<=i || config(i,1)!="1") {
      printf("================================================\n Skip config %02d\n",i); continue; }

    // load configuration
    AliDielectron *jpsi=ConfigJpsi_pd_pp(i,hasMC,runnumber,period);
    if(!jpsi) continue;

    // create unique title
    TString unitit = Form("%s_%s",triggerNames[j],jpsi->GetName());

    // add dielectron to the task and manager
    task->AddDielectron(jpsi);

    printf(" Config %s added\n",jpsi->GetName());

  } //end : loop over configs

  // multiple output connection
  mgr->AddTask(task);

  //create output sub containers
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("pdill_QA", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("pdill_CF", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("pdill_EventStat", TH1D::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  //  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
