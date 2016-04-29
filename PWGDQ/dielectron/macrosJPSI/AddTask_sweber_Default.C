Bool_t isAOD=kFALSE;
Bool_t hasMC=kFALSE;
Int_t iPeriod=-1;
enum { k10b=0, k10c, k10d, k10e, k10f, k10h, k10pp, k11a, k11d, k11h, k12h, k13b, k13c, k13d, k13e, k13f, k15f, k15h };


AliAnalysisTask* AddTask_sweber_Default(
  TString cfg = "ConfigJpsi_sw_pp",
  Bool_t gridconf = kTRUE,
  TString prod = "", 
  Bool_t isMC = kFALSE,
  Bool_t rejectPileup = kTRUE,
  Bool_t usePhysicsSelection = kTRUE,
  TString triggerClass = ""
)
{
  // get the current analysis manager
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sweber_Default", "No analysis manager found.");
    return 0;
  }

  // Do we have an MC handler?
  hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  // AOD input?
  isAOD = mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if(isAOD) hasMC = isMC;

  // Get the current train configuration
  TString trainConfig = gSystem->Getenv("CONFIG_FILE");
  TString list = gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;

  // selected period
  if(      prod.Contains("LHC10b") ) iPeriod = k10b;
  else if( prod.Contains("LHC10c") ) iPeriod = k10c;
  else if( prod.Contains("LHC10d") ) iPeriod = k10d;
  else if( prod.Contains("LHC10e") ) iPeriod = k10e;
  else if( prod.Contains("LHC10f") ) iPeriod = k10f;
  else if( prod.Contains("LHC10h") ) iPeriod = k10h;
  else if( prod.Contains("LHC11a") ) iPeriod = k11a;
  else if( prod.Contains("LHC11d") ) iPeriod = k11d;
  else if( prod.Contains("LHC11h") ) iPeriod = k11h;
  else if( prod.Contains("LHC12h") ) iPeriod = k12h;
  else if( prod.Contains("LHC13b") ) iPeriod = k13b;
  else if( prod.Contains("LHC13c") ) iPeriod = k13c;
  else if( prod.Contains("LHC13d") ) iPeriod = k13d;
  else if( prod.Contains("LHC13e") ) iPeriod = k13e;
  else if( prod.Contains("LHC13f") ) iPeriod = k13f;
  else if( prod.Contains("LHC15f") ) iPeriod = k15f;
  else if( prod.Contains("LHC15h") ) iPeriod = k15h;


  // create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("JpsiDefault");
  if ( usePhysicsSelection ) task->UsePhysicsSelection();


  
  // add special triggers
  switch(iPeriod) {
    case k11d: task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);     break;
    case k11h: task->SetTriggerMask(AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral); break;
    case k12h: task->SetTriggerMask(AliVEvent::kAnyINT); break;                                      
    case k13b: task->SetTriggerMask(AliVEvent::kINT7); break;
    case k13c: task->SetTriggerMask(AliVEvent::kINT7); break;
    case k13d: task->SetTriggerMask(AliVEvent::kAnyINT); break;
    case k13e: task->SetTriggerMask(AliVEvent::kAnyINT); break;
    case k13f: task->SetTriggerMask(AliVEvent::kAnyINT); break;
    case k15f: task->SetTriggerMask(AliVEvent::kINT7); break;
    case k15h: task->SetTriggerMask(AliVEvent::kINT7); break;
  }
  mgr->AddTask(task);
  if(! triggerClass.IsNull() ) task->SetFiredTriggerName(triggerClass.Data() );
  
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
  if (!gROOT->GetListOfGlobalFunctions()->FindObject("ConfigJpsi")){
    gROOT->LoadMacro(configFile.Data());
  }

  ConfigJpsi(task);
  

  //   task->SetTriggerOnV0AND();
  if(rejectPileup) task->SetRejectPileup();

  // create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("sweber_Default_tree",
      TTree::Class(),
      AliAnalysisManager::kExchangeContainer,
      "sweber_Default_default");

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("sweber_Default_QA",
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      "JPSI.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("sweber_Default_CF",
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      "JPSI.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("sweber_Default_EventStat",
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
