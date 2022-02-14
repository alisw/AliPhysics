// Important do not use this FilterTask for final Filtered AOD productions!
// Due to the variability of the used config file, it is to easy to loose track of the used settings!

// ROOT6 modifications
#ifdef __CLING__
  #include <AliAnalysisManager.h>
  #include <AliAODInputHandler.h>
  #include <AliDielectronVarCuts.h>

  // Tell ROOT where to find AliPhysics headers
  R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif
 

AliAnalysisTask *AddTask_jjung_LMEEFilter_PbPb(
  TString cfg="ConfigLMEE_nano_PbPb2015.C",
  Bool_t gridconf=kFALSE,
  TString period="",
  Bool_t storeLS = kFALSE,
  Bool_t hasMC_aod = kFALSE,
  // ULong64_t triggers=AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral+AliVEvent::kEMCEGA+AliVEvent::kEMCEJE
  ULong64_t triggers=AliVEvent::kINT7
  ){

  std::cout << "Start Filtering!!" << std::endl;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLMEEFilter", "No analysis manager found.");
    return 0;
  }
  mgr->SetDebugLevel(10);

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

  TString alienPath("alien:///alice/cern.ch/user/j/jjung/");
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
  /////////// >>>>>>>>> aliroot config
  //else if(!gridconf) configFile=alirootPath.Data();
  /////////// add config to path
  //configFile+="/";
  //configFile+=cfg.Data();

  ///////// >>>>>>>>> local config
  else if(!gridconf) configFile="/data4/jung/nanoAOD/";
  ///////// add config to path
  configFile+="/";
  configFile+=cfg.Data();




  //load dielectron configuration file (only once)
  TString checkconfig=cfg;
  if (gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data())){
     printf("ERROR: config already loaded? %s \n", cfg.Data() );
      return;
  }


  #if defined(__CLING__)
    TMacro conf_die(gSystem->ExpandPathName(configFile.Data())); //ROOT6
    AliDielectron *diEle = reinterpret_cast<AliDielectron *>(conf_die.Exec(Form("0,hasMC,period")));

  #elif defined(__CINT__)
    gROOT->LoadMacro(configFile.Data());
    std::cout << "Config: " << checkconfig << std::endl;
    AliDielectron *diEle=ConfigLMEE_nano_PbPb(0,hasMC,period);

  #endif


  

  if(isAOD) {
    //add options to AliAODHandler to duplicate input event
    AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
    aodHandler->SetCreateNonStandardAOD();
    aodHandler->SetNeedsHeaderReplication();

    aodHandler->SetNeedsTOFHeaderReplication();
    aodHandler->SetNeedsVZEROReplication();
    aodHandler->SetNeedsTrackletsBranchReplication();

    //aodHandler->SetNeedsTracksBranchReplication();
    //aodHandler->SetNeedsCaloClustersBranchReplication();
    //aodHandler->SetNeedsVerticesBranchReplication();
    //aodHandler->SetNeedsCascadesBranchReplication();
    //aodHandler->SetNeedsTrackletsBranchReplication();
    //aodHandler->SetNeedsPMDClustersBranchReplication();
    //aodHandler->SetNeedsJetsBranchReplication();
    //aodHandler->SetNeedsFMDClustersBranchReplication();
    //aodHandler->SetNeedsMCParticlesBranchReplication();
    //aodHandler->SetNeedsDimuonsBranchReplication(); // deactivates several branches
    //    if(hasMC) aodHandler->SetNeedsV0sBranchReplication();
    //if(hasMC) aodHandler->SetNeedsMCParticlesBranchReplication();
    diEle->SetHasMC(hasMC);
  }

  //Create task and add it to the analysis manager
  AliAnalysisTaskDielectronFilter *task=new AliAnalysisTaskDielectronFilter("LMEE_DielectronFilter");
  task->SetTriggerMask(triggers);
  task->SetStoreCaloClusters(kFALSE);
  if (!hasMC) task->UsePhysicsSelection();


  task->SetDielectron(diEle);
  if(storeLS) task->SetStoreLikeSignCandidates(storeLS);
  task->SetCreateNanoAODs(kTRUE);
  //task->SetStoreEventplanes(kTRUE);
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
