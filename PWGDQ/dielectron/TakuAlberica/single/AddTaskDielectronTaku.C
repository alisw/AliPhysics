AliAnalysisTask *AddTaskDielectronTaku(Float_t centrMin, Float_t centrMax, 
				       TString fileName, TString suffixName="", 
				       Bool_t hasMC_aod = kFALSE)
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDielectron", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskDielectron", "This task requires an input event handler");
    return NULL;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(mgr->GetMCtruthEventHandler()!=0x0);
  /*
  TString configFile("./ConfigJpsi2eeDataTaku.C");
  if (hasMC){
    configFile="$ALICE_ROOT/PWG3/dielectron/macros/ConfigJpsi2eeEff.C";
  }
  */
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts(
							       Form("eventCuts_%s",suffixName.Data()),
							       "Vertex Track && |vtxZ|<10 && ncontrib>0"
							       );
  eventCuts->SetCentralityRange(centrMin,centrMax);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);





  //create task and add it to the manager
  //cout<<"AliAnalysisTaskMultiDielectron : "<<configFile.Data()<<endl;
  AliAnalysisTaskMultiDielectronNewTaku *task=new AliAnalysisTaskMultiDielectronNewTaku(
											Form("MultiDie_%s",
											     suffixName.Data()),
											eventCuts
											);



  //load dielectron configuration file
  //cout<<"LoadMacro="<<configFile.Data()<<endl;

  //gROOT->LoadMacro(configFile.Data());
  //cout<<"LoadMacro End="<<configFile.Data()<<endl;
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectronTaku *jpsi=ConfigJpsi2ee(i,isAOD);
    //jpsi->SetPreFilterAllSigns(kTRUE);
    if (isAOD) jpsi->SetHasMC(hasMC_aod);
    if (jpsi) task->AddDielectron(jpsi);
  }

  // add event filter
  task->SetEventFilter(eventCuts);

  // pileup rejection
  task->SetRejectPileup();






  //========= Add tender to the ANALYSIS manager and set default storage =====
  AliTender *tender=new AliTender("AnalysisTender");
  tender->SetCheckEventSelection(kFALSE);
  //tender->SetDefaultCDBStorage("raw://");
  tender->SetDefaultCDBStorage("alien://folder=/alice/data/2011/OCDB");
  //========= Attach TOF supply ======
  AliTOFTenderSupply *TOFtender = new AliTOFTenderSupply("TOFtender");
  TOFtender->SetTOFres(80);
  TOFtender->SetCorrectExpTimes(kFALSE);
  //TOFtender->SetTheorExpTimes(kTRUE);
  ///tender->AddSupply(TOFtender);

  //========= Attach TPC supply ======
  AliTPCTenderSupply *tpcSupply=new AliTPCTenderSupply("TPCtender");
  tpcSupply->SetDebugLevel(2);
  //tpcSupply->SetMip(50.);
  ///tender->AddSupply(tpcSupply);

  ///mgr->AddTask(tender);


  //======== Event plane =============
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB);
  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  ///mgr->AddTask(eventplaneTask);

  ///mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------
  //cout<<"----- "<<mgr->GetCommonFileName()<<" + "<<fileName<<endl;
  TString containerName = fileName+mgr->GetCommonFileName();
  containerName += ":PWG3_dielectron";
    
  //create output container
									  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("jpsi_QA_%s",suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("jpsi_CF_%s",suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("jpsi_EventStat_%s",suffixName.Data()), 
			 TH1D::Class(), 
			 AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("single_tree_%s", suffixName.Data()),
			 //TList::Class(),
			 TTree::Class(),
			 AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  AliAnalysisDataContainer *coutput_ep1 = mgr->CreateContainer("EPStat",
							       TList::Class(), AliAnalysisManager::kOutputContainer,
							       "EventStat_temp.root");


   AliAnalysisDataContainer *coutput_td1 =
     mgr->CreateContainer("tender_event", AliESDEvent::Class(),
			  AliAnalysisManager::kExchangeContainer,"default_tender");
 

  cout<<"containerName.Data = "<<containerName.Data()<<endl;

  //  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput(eventplaneTask,1,coutput_ep1);


  //  mgr->ConnectInput(tender, 0, mgr->GetCommonInputContainer());
  //  mgr->ConnectOutput(tender,1,coutput_td1);


  

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  mgr->ConnectOutput(task, 4, coutput1);



  return task;
}
