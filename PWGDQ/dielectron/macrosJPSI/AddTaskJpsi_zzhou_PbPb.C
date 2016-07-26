Bool_t isAOD=kTRUE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f };


AliAnalysisTask* AddTaskJpsi_zzhou_PbPb(
						TString cfg="ConfigJpsi_zzhou_lowpt",
				    Bool_t gridconf=kTRUE,
						TString prod="",
						Bool_t isMC=kFALSE)
{
  //get the current analysis manager
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskJpsi_zzhou_PbPb", "No analysis manager found.");
    return 0;
  }
  
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //AOD input?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(isAOD) hasMC=isMC;
  
  //Get the current train configuration
  //  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  //TString list=gSystem->Getenv("LIST");
  //if( list.IsNull()) list=prod;
  
  // selected period	   
  if( !prod.CompareTo("LHC10h") ) iPeriod = k10h;
  else if( !prod.CompareTo("LHC11h") ) iPeriod = k11h;
                                       
  //Do 
  // // aod monte carlo
  // if( list.Contains("LHC11a10") ||
  //     list.Contains("LHC11b10") ||
  //     list.Contains("LHC12a17") ||
  //     list.Contains("fix")
  //     ) hasMC=kTRUE;
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("Jpsi_zzhouPbPb");

  if (!hasMC) task->UsePhysicsSelection();
  
  // add special triggers
  switch(iPeriod) {
  case k11h: task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral); break;
  }
  mgr->AddTask(task);
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
    //  task->SetTriggerOnV0AND();
  task->SetEventFilter(eventCuts);
  task->SetRejectPileup();
		
  gROOT->LoadMacro( cfg.Data() );
	
	
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_zzhou_lowpt(i);
    if (!jpsi) continue;
    jpsi->SetHasMC(hasMC);
    task->AddDielectron(jpsi);
  }
  
  //   task->SetTriggerOnV0AND();
  //   if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("jpsi_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "jpsi_Default_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "AnalysisResults.root");
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "AnalysisResults.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "AnalysisResults.root");
  
	mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 0, coutput1 );
	mgr->ConnectOutput(task, 1, cOutputHist1);
	mgr->ConnectOutput(task, 2, cOutputHist2);
	mgr->ConnectOutput(task, 3, cOutputHist3);
	
	return task;
}
