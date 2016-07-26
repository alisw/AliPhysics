AliAnalysisTask* AddTask_sweber_Default(
  TString cfg = "ConfigJpsi_sw_pp_13tev",
  Bool_t gridconf = kTRUE,
  TString suffix = ""
)
{
  // get the current analysis manager
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sweber_Default", "No analysis manager found.");
    return 0;
  }

  
  // create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron( Form("JpsiDefault%s",suffix.Data()) );
  mgr->AddTask(task);

 
  // set config file name
  TString configFile("");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsi_sw_pp";

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
  //}

  ConfigJpsi(task);

  // create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer( Form("sweber_Default_tree%s",suffix.Data()),
      TTree::Class(),
      AliAnalysisManager::kExchangeContainer,
      Form("sweber_Default_default%s",suffix.Data()));

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("sweber_Default_QA%s",suffix.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("JPSI%s.root",suffix.Data()));

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("sweber_Default_CF%s",suffix.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("JPSI%s.root",suffix.Data()));

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("sweber_Default_EventStat%s",suffix.Data()),
      TH1D::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("JPSI%s.root",suffix.Data()));

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;
}
