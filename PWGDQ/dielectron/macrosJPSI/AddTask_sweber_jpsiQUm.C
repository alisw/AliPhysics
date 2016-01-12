AliAnalysisTask *AddTask_sweber_jpsiQUm(
  TString cfg = "ConfigJpsiTriggerQU",
  Bool_t gridconf = kTRUE,
  TString prod = "", 
  Bool_t isMC = kFALSE,
  Bool_t rejectPileup = kTRUE){
	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sweber_jpsiQUm", "No analysis manager found.");
    return 0;
  }


  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron("sweberTaskTRD");
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if (!hasMC) task->UsePhysicsSelection(  );
  task->SetTriggerMask(AliVEvent::kTRD);
 // task->SetTriggerMask(AliVEvent::kINT7);
 // if (!hasMC) task->SelectCollisionCandidates(AliVEvent::kTRD);

 
	task->SetRequireTRDTrigger( kTRUE);
	task->SetTRDTriggerClass(AliDielectronEventCuts::kQU);
	task->SetRequireMatchedTrack(1);
 
 /**
 *
 * train config
 *
 **/
 
  
  // set config file name
  TString configFile("");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsiTriggerQU";

  // the different paths
  TString gsiPath("$TRAIN_ROOT/sweber_chic");
  TString alienPath("alien:///alice/cern.ch/user/s/sgweber/macrosJPSI");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI");

  // >>> gsi config
  if (!trainRoot.IsNull()){
    configFile=gsiPath.Data();
  }
  // >>> aliroot config
  else if(!gridconf && trainRoot.IsNull()){
    configFile=alirootPath.Data();
  }
  // >>> alien config
  else{
    if(!gSystem->Exec(Form("alien_cp %s/%s.C .",alienPath.Data(),cfg.Data()))) {
          configFile=gSystem->pwd();
    }
    else {
      printf("ERROR: couldn't copy file %s/%s.C from grid \n", alienPath.Data(),cfg.Data() );
      return;
    }
  }
  // add config to path
  configFile+="/";
  configFile+=cfg.Data();
  configFile+=".C";
  
  // load dielectron configuration file (only once)
 // if (!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigJpsi")){
    gROOT->LoadMacro(configFile.Data());
 // }

  ConfigJpsi(task);
  

  //   task->SetTriggerOnV0AND();
  if(rejectPileup) task->SetRejectPileup();




	mgr->AddTask(task);
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("sweber_treeQUm",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "sweber_defaultQUm");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsiQA_QUm",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsiCF_QUm",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsiEventStat_QUm",
      TH1D::Class(),
      AliAnalysisManager::kOutputContainer,
      "JPSI.root");

  AliAnalysisDataContainer *cOutputHist4 =
    mgr->CreateContainer("jpsiEventStatTRDtrigger_QUm",
      TH1D::Class(),
      AliAnalysisManager::kOutputContainer,
      "JPSI.root");

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  mgr->ConnectOutput(task, 4, cOutputHist4);

  return task;
}
