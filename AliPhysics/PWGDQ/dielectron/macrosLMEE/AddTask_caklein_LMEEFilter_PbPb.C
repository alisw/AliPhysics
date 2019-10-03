// Important do not use this FilterTask for final Filtered AOD productions!
// Due to the variability of the used config file, it is to easy to loose track of the used settings!

AliAnalysisTask *AddTask_caklein_LMEEFilter_PbPb(
  TString cfg="ConfigLMEE_nano_PbPb2015.C",
  Bool_t gridconf=kFALSE,
  TString period="",
  Bool_t storeLS = kFALSE,
  Bool_t hasMC_aod = kFALSE,
  // ULong64_t triggers=AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral+AliVEvent::kEMCEGA+AliVEvent::kEMCEJE
  ULong64_t triggers=AliVEvent::kINT7
  ){

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


  TString configFile("");
  printf("pwd:        %s \n",gSystem->pwd());
  if(cfg.IsNull()) cfg="ConfigLMEE_nano_PbPb2015.C";

  TString alienPath("alien:///alice/cern.ch/user/c/cklein/PWGDQ/dielectron/macrosLMEE/");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");

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

  //load dielectron configuration file (only once)
  TString checkconfig="ConfigLMEE_nano_PbPb2015";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  AliDielectron *diEle=ConfigLMEE_nano_PbPb2015(0,hasMC,period);

  if(isAOD) {
    //add options to AliAODHandler to duplicate input event
    AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
    aodHandler->SetCreateNonStandardAOD();
    aodHandler->SetNeedsHeaderReplication();

   if(!period.Contains("LHC10h")) aodHandler->SetNeedsTOFHeaderReplication();
    aodHandler->SetNeedsVZEROReplication();
    /*aodHandler->SetNeedsTracksBranchReplication();
    aodHandler->SetNeedsCaloClustersBranchReplication();
    aodHandler->SetNeedsVerticesBranchReplication();
    aodHandler->SetNeedsCascadesBranchReplication();
    aodHandler->SetNeedsTrackletsBranchReplication();
    aodHandler->SetNeedsPMDClustersBranchReplication();
    aodHandler->SetNeedsJetsBranchReplication();
    aodHandler->SetNeedsFMDClustersBranchReplication();
    //aodHandler->SetNeedsMCParticlesBranchReplication();
    aodHandler->SetNeedsDimuonsBranchReplication();*/ // deactivates several branches
    //    if(hasMC) aodHandler->SetNeedsV0sBranchReplication();
    if(hasMC) aodHandler->SetNeedsMCParticlesBranchReplication();
    diEle->SetHasMC(hasMC);
  }

  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("LMEE_DielectronFilter");
  task->SetTriggerMask(triggers);
  if (!hasMC) task->UsePhysicsSelection();


  task->SetDielectron(diEle);
  if(storeLS) task->SetStoreLikeSignCandidates(storeLS);
  task->SetCreateNanoAODs(kTRUE);
  task->SetStoreEventplanes(kTRUE);
  task->SetStoreEventsWithSingleTracks(kTRUE); 
  // task->SetStoreHeader(kTRUE);
  mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------


  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGDQ_dielectronFilter";

  //create output container

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("LMEE_FilterQA",
                         THashList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("LMEE_FilterEventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());


  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);

  return task;
}
