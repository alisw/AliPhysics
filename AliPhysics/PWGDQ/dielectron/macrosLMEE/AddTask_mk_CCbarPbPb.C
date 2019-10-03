AliAnalysisTask *AddTask_mk_CCbarPbPb(TString prod="",ULong64_t triggers=AliVEvent::kEMCEGA  | AliVEvent::kEMCEJE){
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_mk_CCbarPbPb", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_mk_CCbarPbPb", "This task requires an input event handler");
    return NULL;
  }

     
    //Do we have an MC handler?
  Bool_t hasMC = kFALSE;
  TString list = gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;
  
  
TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/Config_mk_CCbarPbPb.C");
//TString configFile("Config_mk_CCbarPbPb.C");
 
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //create task and add it to the manager
   AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDie");

     //load dielectron configuration file
TString checkconfig="$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/Config_mk_CCbarPbPb";
//  TString checkconfig="ConfigCCbar_mk_pp";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi; 
    jpsi=Config_mk_CCbarPbPb(i,isAOD);
    if (jpsi) task->AddDielectron(jpsi);
  }

  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetCentralityRange(0.0,90.0);
  // add event filter
  task->SetEventFilter(eventCuts);

  // pileup rejection
  task->SetTriggerMask(AliVEvent::kEMCEGA);
  task->UsePhysicsSelection();
   
  mgr->AddTask(task);
  
  //----------------------
  //create data containers
  //----------------------

  //create output container
  TString containerName = "LMEEoutput.root";
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("mkoehler_QA_EMC",
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("mkoehler_CF_EMC",
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("mkoehler_EventStat_EMC",
			 TH1D::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());
  
 
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
