

AliAnalysisTask *AddTask_mkoehler_CCbarMC(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_mkoehler_CCbarMC", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_mkoehler_CCbarMC", "This task requires an input event handler");
    return NULL;
  }

TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/Config_mkoehler_CCbarMC.C");
//TString configFile("Config_mkoehler_CCbarMC.C");


  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron("MultiDie");
  if (!hasMC) return ;



  TString checkconfig="$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/Config_mkoehler_CCbarMC";
//  TString checkconfig="ConfigCCbar_mkoehler_CCbarMC";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  task->SetEventFilter(eventCuts);

//  task->SetRejectPileup();

  mgr->AddTask(task);
  
  //load dielectron configuration file
//  gROOT->LoadMacro("Config_mkoehler_CCbarMC.C");
//gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/Config_mkoehler_CCbarMC.C");

  //add dielectron analysis with different cuts to the task

  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *diel_lowMC = Config_mkoehler_CCbarMC(i,isAOD);
    task->AddDielectron(diel_lowMC);

  }
  

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmassMC",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_diel_lowmassMC",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_diel_lowmassMC",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");
  

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("mkoehler_lowmassMC_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
