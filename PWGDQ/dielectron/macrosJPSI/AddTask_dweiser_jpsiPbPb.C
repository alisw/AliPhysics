AliAnalysisTask *AddTask_dweiser_jpsiPbPb( 
   TString cfg = "ConfigJpsi_dw_Conversions",
   Bool_t gridconf = kFALSE
)
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_dweiser_jpsiPbPb", "No analysis manager found.");
    return 0;
  }

  //TString configFile("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_dw_Conversions.C");

  // TString configFile("");
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("JpsiDefault");
  mgr->AddTask(task);

  TString configFile("");

  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsi_dw_Conversions";



  // the different paths                                                                                                                       
  TString gsiPath("$TRAIN_ROOT/dweiser_ff/jpsiPbPb/src");
  TString alienPath("alien:///alice/cern.ch/user/d/dweiser/macrosJPSI");
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



  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  // for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
  //  AliDielectron *jpsi=ConfigDefault(i, task);
  //  if (!jpsi) continue;
  //  task->AddDielectron(jpsi);
  // }

  ConfigDefault(task);

  //Event cuts applied in Config
  /* 
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //  eventCuts->SetCentralityRange(0.0,100.0); //does sth?
  task->SetEventFilter(eventCuts);
 
  //   task->SetTriggerOnV0AND();
   if ( trainConfig=="pp" ) task->SetRejectPileup();
  */
  //create output container


  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("jpsi_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "out_dweiser_jpsiPbPb_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "out_dweiser_jpsiPbPb.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "out_dweiser_jpsiPbPb.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "out_dweiser_jpsiPbPb.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
