AliAnalysisTask *AddTask_jbook_v2(TString prod="",
				    Bool_t gridconf=kFALSE,
				    ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB) {

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jbook_v2", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC = kFALSE;
  TString list = gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;

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
  if (!trainRoot.IsNull())
    configFile="$TRAIN_ROOT/jbook_jpsi/ConfigJpsi_jb_PbPb.C";   // gsi config
  else if(!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/j/jbook/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_jb_PbPb.C ."))
    configFile=Form("%s/ConfigJpsi_jb_PbPb.C",gSystem->pwd());                        // alien config

  if(!gridconf)
    configFile="$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_jb_PbPb.C"; // aliroot config

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task;

  // trigger selection
  ULong64_t triggerSets[]={AliVEvent::kCentral , AliVEvent::kSemiCentral , AliVEvent::kMB,
			   AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB};
  const char* triggerNames[]={"Central","SemiCentral","MB","MB+Cent+SemiCent"};

  // find out the configured triggers
  Int_t j=0;
  for(j=0; j<4; j++) {
    if(triggers!=triggerSets[j]) continue;
    else break;
  }

  // print task configuration
  printf("production: %s MC: %d \n",  list.Data(),hasMC);
  printf("triggers:   %s \n",         triggerNames[j]  );
  printf("config:     %s Grid: %d \n",configFile.Data(),gridconf);

  task = new AliAnalysisTaskMultiDielectron((Form("MultiDieJB_%s",triggerNames[j])));
  task->SetTriggerMask(triggers);
  //task->SetTriggerMask(AliVEvent::kMB);

  if (!hasMC) task->UsePhysicsSelection();
  mgr->AddTask(task);

  //load dielectron configuration file
  TString checkconfig="ConfigJpsi_jb_PbPb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  //add dielectron analysis with different cuts to the task
  for (Int_t i=12; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_jb_PbPb(i,list,triggers);
    if (jpsi ) task->AddDielectron(jpsi);
    if (jpsi ) printf(" %s added\n",jpsi->GetName());
  }

  //create output container
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("jbook_v2_%s",triggerNames[j]),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("jbook_v2_%s.root",triggerNames[j]));

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("jbook_v2_CF_%s",triggerNames[j]),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("jbook_v2_%s.root",triggerNames[j]));

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("jbook_v2_EventStat_%s",triggerNames[j]),
			 TH1D::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("jbook_v2_%s.root",triggerNames[j]));

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  //  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
