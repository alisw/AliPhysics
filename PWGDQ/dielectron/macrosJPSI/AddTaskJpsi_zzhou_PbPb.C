Bool_t isAOD=kTRUE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10h, k11h };


AliAnalysisTask* AddTaskJpsi_zzhou_PbPb(
					   TString cfg="ConfigJpsi_zzhou_lowpt.C",
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

  //set config file name
  TString configFile("");
  printf("pwd:        %s \n",gSystem->pwd());
  if(cfg.IsNull()) cfg="ConfigJpsi_zzhou_lowpt.C";

  // the different paths
  TString alienPath("alien:///alice/cern.ch/user/z/zzhou/PWGDQ/dielectron/macrosJPSI/");
  TString alirootPath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosJPSI/");

  // #####  alien config  #####
  if(gridconf) {
    if(!gSystem->Exec(Form("alien_cp %s/%s .",alienPath.Data(),cfg.Data()))) {
      gSystem->Exec(Form("ls -l %s",gSystem->pwd()));
      configFile=gSystem->pwd();
    }
    else {
      printf("ERROR: couldn't copy file %s/%s from grid \n", alienPath.Data(),cfg.Data() );
      return;
    }
  }
  
  // ##### aliroot config #####
  else if(!gridconf) configFile=alirootPath.Data();
  // ##### add config to path #####
  configFile+="/";
  configFile+=cfg.Data();  

  // print overall configuration
  printf("config:     %s Grid: %d \n",configFile.Data(),gridconf);
  printf("#######------------overall configuration----------------#########\n");
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("Jpsi_zzhouPbPb");

  if (!hasMC) task->UsePhysicsSelection();
  
  // add special triggers
  switch(iPeriod) {
  case k11h: task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral); break;
  }


  // event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
    //  task->SetTriggerOnV0AND();
  eventCuts->SetCentralityRange(40.,90.);
  task->SetEventFilter(eventCuts);
  task->SetRejectPileup();
		
  // use histogram to correct PID need etaCorr_2010.root from alien
  //TFile *file =TFile::Open("alien:///alice/cern.ch/user/z/zzhou/PWGDQ/dielectron/files/etaCorr_2010.root");

  //load dielectron configuration file (only once)
  TString checkconfig="ConfigJpsi_zzhou_lowpt";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  //define default output container
  TString containerName = "AnalysisResults.root";
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi_zzhou_lowpt(i);
    if (!jpsi) continue;
    jpsi->SetHasMC(hasMC);
    task->AddDielectron(jpsi);
  } //end : loop over configs
  
  mgr->AddTask(task);
  
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
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  
	mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 0, coutput1 );
	mgr->ConnectOutput(task, 1, cOutputHist1);
	mgr->ConnectOutput(task, 2, cOutputHist2);
	mgr->ConnectOutput(task, 3, cOutputHist3);
	
	return task;
}
