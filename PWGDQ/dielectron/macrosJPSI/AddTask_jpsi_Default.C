Bool_t isAOD=kFALSE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f };


AliAnalysisTask* AddTask_jpsi_Default(
						TString cfg="ConfigJpsi_raul_qa",
				    Bool_t gridconf=kTRUE,
						TString prod="",
						Bool_t isMC=kFALSE)
{
  //get the current analysis manager
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_Default", "No analysis manager found.");
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
  
  //printf("list************* %s",list.Data());
  // selected period
	
  if(      !prod.CompareTo("LHC10b") ) iPeriod = k10b;
  else if( !prod.CompareTo("LHC10c") ) iPeriod = k10c;
  else if( !prod.CompareTo("LHC10d") ) iPeriod = k10d;
  else if( !prod.CompareTo("LHC10e") ) iPeriod = k10e;
  else if( !prod.CompareTo("LHC10f") ) iPeriod = k10f;
  else if( !prod.CompareTo("LHC10h") ) iPeriod = k10h;
  else if( !prod.CompareTo("LHC11a") ) iPeriod = k11a;
  else if( !prod.CompareTo("LHC11d") ) iPeriod = k11d;
  else if( !prod.CompareTo("LHC11h") ) iPeriod = k11h;
  else if( !prod.CompareTo("LHC12h") ) iPeriod = k12h;
  else if( !prod.CompareTo("LHC13b") ) iPeriod = k13b;
  else if( !prod.CompareTo("LHC13c") ) iPeriod = k13c;
  else if( !prod.CompareTo("LHC12d") ) iPeriod = k13d;
  else if( !prod.CompareTo("LHC12e") ) iPeriod = k13e;
  else if( !prod.CompareTo("LHC12f") ) iPeriod = k13f;
  

  //
  //////                                                                                                                                                                                                                                             
  //Do 
  // // aod monte carlo
  // if( list.Contains("LHC11a10") ||
  //     list.Contains("LHC11b10") ||
  //     list.Contains("LHC12a17") ||
  //     list.Contains("fix")
  //     ) hasMC=kTRUE;
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("JpsiDefault");
  //  task->SetBeamEnergy(1380.); // not neeeded since we are not looking at helicity and Collins-Soper coordinates
  if (!hasMC) task->UsePhysicsSelection();
  
  // add special triggers
  switch(iPeriod) {

  case k10b: task->SetTriggerMask(AliVEvent::kMB); break;
  case k10c: task->SetTriggerMask(AliVEvent::kMB); break;
  case k10d: task->SetTriggerMask(AliVEvent::kMB); break;
  case k10e: task->SetTriggerMask(AliVEvent::kMB); break;
  case k11d: task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);     break;
  case k11h: task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral); break;
  case k12h: task->SetTriggerMask(AliVEvent::kAnyINT); break;                                      
  case k13b: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k13c: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k13d: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  case k13e: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  case k13f: task->SetTriggerMask(AliVEvent::kAnyINT); break;
  case k15f: task->SetTriggerMask(AliVEvent::kINT7); break;
  case k15n: task->SetTriggerMask(AliVEvent::kINT7); break;


  }
  mgr->AddTask(task);
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.); 
  //  task->SetTriggerOnV0AND();
  task->SetEventFilter(eventCuts);
  // task->SetRejectPileup();
	
	
		  //set config file name
  TString configFile("");
  printf("%s \n",gSystem->pwd());
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if(cfg.IsNull()) cfg="ConfigJpsi_raul_qa";

  // the different paths
  TString gsiPath("$TRAIN_ROOT/rjimenez_jpsi");
  TString alienPath("alien:///alice/cern.ch/user/r/raul/macrosJPSI");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI");

  ////////// >>>>>>>>>> gsi config
  if (!trainRoot.IsNull()){
		configFile=gsiPath.Data();
	}
  ///////// >>>>>>>>> aliroot config
	else if(!gridconf && trainRoot.IsNull()){
		configFile=alirootPath.Data();
	}
  ////////// >>>>>>>>>> alien config
  else{
		if(!gSystem->Exec(Form("alien_cp %s/%s.C .",alienPath.Data(),cfg.Data()))) {
    gSystem->Exec(Form("ls -l %s",gSystem->pwd()));
    configFile=gSystem->pwd();
		}
		else {
			printf("ERROR: couldn't copy file %s/%s.C from grid \n", alienPath.Data(),cfg.Data() );
			return;
		}
	}
  
  ///////// add config to path
  configFile+="/";
  configFile+=cfg.Data();
  configFile+=".C";
	
	  //load dielectron configuration file (only once)
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(cfg.Data())){
    gROOT->LoadMacro(configFile.Data());
	}
	

	
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_raul_qa(i);
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
                         "JPSI.root");
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "JPSI.root");
  
	mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 0, coutput1 );
	mgr->ConnectOutput(task, 1, cOutputHist1);
	mgr->ConnectOutput(task, 2, cOutputHist2);
	mgr->ConnectOutput(task, 3, cOutputHist3);
	
	return task;
}
