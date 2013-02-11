AliAnalysisTask *AddTask_mfigueredo_JPsi(TString prod="",ULong64_t triggers=AliVEvent::kEMCEGA  | AliVEvent::kEMCEJE){
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_mfigueredo_JPsi", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_mfigueredo_JPsi", "This task requires an input event handler");
    return NULL;
  }

     
    //Do we have an MC handler?
  Bool_t hasMC = kFALSE;
  TString list = gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;
  
  TString configFile("");
  printf("%s \n",gSystem->pwd());
  
  configFile="ConfigJpsi_mf_PbPb.C";
 
  if(!gSystem->Exec("alien_cp alien:///alice/cern.ch/user/m/mfiguere/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mf_PbPb.C ."))
    configFile=Form("%s/ConfigJpsi_mf_PbPb.C",gSystem->pwd());                        // alien config                                                                                            
  else
    configFile="$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigJpsi_mf_PbPb.C"; // aliroot config                                               
  
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //create task and add it to the manager
//   AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDie");
 
  // trigger selection
  ULong64_t triggerSets[]={AliVEvent::kEMCEGA ,AliVEvent::kEMCEJE};
  const char* triggerNames[]={"EMCEGA","EMCEJE"};

  // find out the configured triggers
  Int_t j=0;
  for(j=0; j<2; j++) {
    if(triggers!=triggerSets[j]) continue;
    else break;
  }
  
    // print task configuration
  printf("production: %s MC: %d \n",  list.Data(),hasMC);
  printf("triggers:   %s \n",         triggerNames[j]  );

  task = new AliAnalysisTaskMultiDielectron((Form("MultiDieData_%s",triggerNames[j])));
   
  //load dielectron configuration file
  TString checkconfig="ConfigJpsi_mf_pp";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi; 
    jpsi=ConfigJpsi_mf_PbPb(i,isAOD);
    if (jpsi) task->AddDielectron(jpsi);
    if (jpsi ) printf(" %s added\n",jpsi->GetName());
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
  task->SetTriggerMask(triggers);
  task->UsePhysicsSelection();
   
  mgr->AddTask(task);
  
  //----------------------
  //create data containers
  //----------------------

  //create output container
  TString containerName = "JPSI.root";
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("mfigueredo_QA_%s",triggerNames[j]),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("mfigueredo_CF_%s",triggerNames[j]),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("mfigueredo_EventStat_%s",triggerNames[j]),
			 TH1D::Class(),
			 AliAnalysisManager::kOutputContainer,
			 containerName.Data());
  
 
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
