AliAnalysisTask *AddTaskLMEEFilter_PbPb(TString cfg="ConfigLMEE_nano_PbPb.C",
				      Bool_t gridconf=kFALSE,
				      ULong64_t triggers=AliVEvent::kINT7,
				      TString period="",
				      Bool_t useTrackCuts = kTRUE,
				      Bool_t storeLS = kTRUE,
				      Bool_t hasMC_aod = kFALSE){

  // This AddTask macro is a copy of the JPsi Nano Filteringtask:
  // (PWGDQ/dielectron/macrosJPSI/AddTaskJPSIFilter_pp.C
  // and adapted to LMEE analyses

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLMEEFilter", "No analysis manager found.");
    return 0;
  }
  
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Warning("AddTaskLMEEFilter","No AOD output handler available. Not adding the task!");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0)||hasMC_aod;
  
  //Do we run on AOD?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //Allow merging of the filtered aods on grid trains
  if(mgr->GetGridHandler()) {
    printf(" SET MERGE FILTERED AODs \n");
    //mgr->GetGridHandler()->SetMergeAOD(kTRUE);
  }

  //set config file name
  TString configFile("");
  printf("%s \n",gSystem->pwd());
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigLMEE_nano_PbPb.C";

  // the different paths
  TString alienPath("alien:///alice/cern.ch/user/m/miweber/PWGDQ/dielectron/macrosLMEE");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");

  ////////// >>>>>>>>>> alien config
  if(gridconf){
    if(!gSystem->Exec(Form("alien_cp %s/%s .",alienPath.Data(),cfg.Data()))){
      gSystem->Exec(Form("ls -l %s",gSystem->pwd()));
      configFile=gSystem->pwd();
    }
    else {
      printf("ERROR: couldn't copy file %s/%s from grid \n", alienPath.Data(),cfg.Data() );
      return;
    }
  }
  else{
    configFile=alirootPath.Data();
  }
  ///////// >>>>>>>>> aliroot config

  ///////// add config to path
  configFile+="/";
  configFile+=cfg.Data();

  //load dielectron configuration file (only once)
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(cfg.Data()))
    gROOT->LoadMacro(configFile.Data());

  AliDielectron *lmee=ConfigLMEE_nano_pp(0,hasMC,period,useTrackCuts);
  
  if(isAOD) {
    //add options to AliAODHandler to duplicate input event
    AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
    aodHandler->SetCreateNonStandardAOD();
    aodHandler->SetNeedsHeaderReplication();
    if(!period.Contains("LHC10h")) aodHandler->SetNeedsTOFHeaderReplication();
    aodHandler->SetNeedsVZEROReplication();
    if(hasMC) aodHandler->SetNeedsMCParticlesBranchReplication();
    lmee->SetHasMC(hasMC);
  }
  
  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("lmee_DielectronFilter");
  task->SetTriggerMask(triggers);
   if (!hasMC) task->UsePhysicsSelection();

  task->SetDielectron(lmee);
  if(storeLS) task->SetStoreLikeSignCandidates(storeLS);
  task->SetCreateNanoAODs(kTRUE);
  task->SetStoreEventsWithSingleTracks(kTRUE);
  mgr->AddTask(task);


  //----------------------
  //create data containers
  //----------------------
    
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGDQ_dielectronFilter";
 
  //create output container
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("lmee_FilterQA",
                         THashList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("lmee_FilterEventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  
  return task;
}
