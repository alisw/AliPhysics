Bool_t isAOD=kTRUE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f,k15f,k15o,k17q,k17p};

AliAnalysisTask* AddTask_xbai_pp_QA(TString cfg="ConfigJpsi_xbai_pp_QA",
				   Bool_t gridconf=kFALSE,
				   TString prod="LHC10h",
				   Bool_t isMC=kFALSE){
  //get the current analysis manager
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_xbai_pp_QA", "No analysis manager found.");
    return 0;
  }
  
  //Do we have an MC handler?
  hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //AOD input?
  isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(isAOD) hasMC=isMC;
  
  
  // selected period	   
  if( !prod.CompareTo("LHC10h") ) iPeriod = k10h;
  else if( !prod.CompareTo("LHC11h") ) iPeriod = k11h;
  else if( !prod.CompareTo("LHC15o") ) iPeriod = k15o;
  else if( !prod.CompareTo("LHC17p") ) iPeriod = k17p;
  else if( !prod.CompareTo("LHC17q") ) iPeriod = k17q;
                               
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("Jpsi_pp");

  if (!hasMC) task->UsePhysicsSelection();

 //Get the current train configuration
  //TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  TString list=gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;
//Do 
  // // aod monte carlo
  if( list.Contains("LHC11a10") ||
      list.Contains("LHC11b10") ||
      list.Contains("LHC12a17") ||
      list.Contains("fix")
      ) hasMC=kTRUE;
  
  
  // add special triggers
  switch(iPeriod) {
  case k10h: task->SetTriggerMask(AliVEvent::kMB);
  case k11h: task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral); break;
  case k15o: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k17p: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k17q: task->SetTriggerMask(AliVEvent::kINT7); break;
  }
  mgr->AddTask(task);

  
// set config file name
  TString configFile("");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsi_xbai_pp_QA";
  // the different paths
  TString gsiPath("$TRAIN_ROOT/rjimenez_jpsi");
  TString alienPath("alien:///alice/cern.ch/user/x/xbai/pp2017/check_nsigmaE");
  //TString alirootPath("/u/xbai/work/pp2017/check_nsigmaE/Train_test");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI");

  //printf("test locall XXX");
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
  printf("this is the file %s xxx",configFile.Data());
    // load dielectron configuration file (only once)
  if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config")){
    gROOT->LoadMacro(configFile.Data());
  }
  Config(task);

  /*
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
    //  task->SetTriggerOnV0AND();
  task->SetEventFilter(eventCuts);
  task->SetRejectPileup();
		
  */	
       
  
  //   task->SetTriggerOnV0AND();
  //   if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container`
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("rjimenez_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "rjimenez_Default_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("rjimenez_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSITree.root");
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("rjimenez_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSITree.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("rjimenez_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSITree.root");
  
	mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 0, coutput1 );
	mgr->ConnectOutput(task, 1, cOutputHist1);
	mgr->ConnectOutput(task, 2, cOutputHist2);
	mgr->ConnectOutput(task, 3, cOutputHist3);
	
	return task;
}
