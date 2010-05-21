//
// Macro to run performance QA train
// locally. The TPC performance task is attached.
// 
//
//13.10.2009 -  J.Otwinowski@gsi.de
//
//

/*
Quick Start:

1. Start train macro (real data from GSI::SE). Source your own Alien environment.

source /u/jacek/alien/set_alien_env.sh
aliroot -b -q 'RunPerformanceTrain.C("AliESDs.root",2,kFALSE,kTRUE,kTRUE)'

3. Start train macro (real data from lustre)

aliroot -b -q 'RunPerformanceTrain.C("AliESDs.root",2,kFALSE,kTRUE,kFALSE)'

*/

//_____________________________________________________________________________
void RunPerformanceTrain(Char_t *file="esd.root", Int_t runNumber = 2, const char* triggerClass ="CINT1B-ABCE-NOPF-ALL", Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, Bool_t bGrid=kFALSE)
{
  //
  // Grid settings
  // use GSI::SE
  if(bGrid) {
    gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE");
    TGrid * alien = TGrid::Connect("alien://",0,0,"t");
    gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE");
  }

  //
  // Train Configuration
  //
  Int_t       iPWG1perfTPC       = 1;      // Test TPC performance
  Int_t       iPWG1perfTRD       = 0;      // Test TRD performance
  Int_t       iPWG1perfITS       = 0;      // Test ITS performance
  Int_t       iPWG1perfCalo      = 0;      // Test Calo performance
  Int_t       iPWG1perfMuonTrig  = 0;      // Test Muon Trigger performance
  Int_t       iPWG1perfMuonEff   = 0;      // Test Muon Efficiency performance
  Int_t       iPWG1perfTOF       = 0;      // Test TOF-TPC matching performance
  Int_t       iPWG1perfPrimVertex = 0;     // Test Primary Vertex performance
  Int_t       iPWG1v0QA          = 0;      // V0 algorithm QA task

  //
  // Load Libraries
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");

  gSystem->Load("libTPCcalib.so");
  gSystem->Load("libPWG1");

  gSystem->Load("libPHOSUtils");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPWG4PartCorrBase");
  gSystem->Load("libPWG4PartCorrDep");

  gSystem->Load("libPWG3muon.so"); // The class is here

  //
  // OCDB Configuration 
  //
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local:///lustre/alice/alien/alice/data/2009/OCDB");
  //cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //cdbManager->SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  cdbManager->SetRun(runNumber);
  //cdbManager->SetCacheFlag(kFALSE);
  // initialize magnetic field from the GRP manager.
  //if(magField==0) TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 0., 0., AliMagF::k2kG));
  //if(magField==1) TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k2kG));
  //if(magField==2) TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

  /*
  AliGRPManager grpMan;
  grpMan.ReadGRPEntry();
  grpMan.SetMagField();
  AliRunInfo *runInfo = grpMan.GetRunInfo();
 
  //
  // Load geometry
  //
  */
  AliGeomManager::LoadGeometry();

  //
  // Swtich off all AliInfo (too much output!)
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //
  // Create input ESD chain
  //
  /*
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(list,nFiles,fistFile);
  if(!chain) { 
    Error("RunPerformanceTrain","ESD chain not created!");
    return;
  }
  */
  TChain  *chain = new TChain("esdTree");
  if(!chain) { 
    Error("RunPerformanceTrain","ESD chain not created!");
    return;
  }
  chain->Add(file);
  chain->Lookup();

  //
  // Create analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager;
  if(!mgr) { 
    Error("RunPerformanceTrain","AliAnalysisManager not set!");
    return;
  }

  //
  // Set ESD input handler
  //
  AliESDInputHandler* esdH = new AliESDInputHandler;
  if(!esdH) { 
    Error("RunPerformanceTrain","AliESDInputHandler not created!");
    return;
  }
  if(bUseESDfriend) esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);

  /*
  //
  // Set RecPoints and ESD input handler
  //
  AliESDInputHandlerRP* rpH = new AliESDInputHandlerRP;
  if(!rpH) { 
    Error("RunPerformanceTrain","AliESDInputHandlerRP not created!");
    return;
  }
  if(bUseESDfriend) rpH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(rpH);
  */

  //
  // Set MC input handler
  //
  if(bUseMCInfo) {
    AliMCEventHandler* mcH = new AliMCEventHandler;
    if(!mcH) { 
      Error("RunPerformanceTrain","AliMCEventHandler not created!");
      return;
    }
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH);
  }


  //
  // Add task to AliAnalysisManager
  //

  //
  // TPC performance
  //
  if(iPWG1perfTPC) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPC.C");
    AliPerformanceTask *tpcQA = AddTaskPerformanceTPC(bUseMCInfo,bUseESDfriend,triggerClass);
    if(!tpcQA) { 
      Error("RunPerformanceTrain","AliPerformanceTask not created!");
      return;
    }
  } 
  else {
    Printf("RunPerformanceTrain: TPC Performance - EXCLUDED!");
  }
  //
  // TRD perormance
  //
  if(iPWG1perfTRD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTrainPerformanceTRD.C");
    if(!AddTrainPerformanceTRD(bUseMCInfo,bUseESDfriend)) { 
      Error("RunPerformanceTrain","TrainPerformanceTRD not created!");
      return;
    }
  } else {
    Printf("RunPerformanceTrain: TRD TrainPerformanceTRD - EXCLUDED!");
  }
  //
  // ITS performance
  //
  if(iPWG1perfITS) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceITS.C");
    AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(bUseMCInfo);
    if(!itsQA) { 
      Error("RunPerformanceTrain","AliAnalysisTaskITSTrackingCheck not created!");
      return;
    }
  } 
  else {
    Printf("RunPerformanceTrain: AliAnalysisTaskITSTrackingCheck - EXCLUDED!");
  }
  //
  // Calorimeter Performance
  //
  if(iPWG1perfCalo) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskCalorimeterQA.C");
    AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD",bUseMCInfo,kFALSE);
    if(!taskCaloQA) { 
      Error("RunPerformanceTrain","AliAnalysisTaskParticleCorrelation not created!");
      return;
    }
    mgr->AddTask(taskCaloQA);
  } 
  else {
    Printf("RunPerformanceTrain: AliAnalysisTaskParticleCorrelation - EXCLUDED!");
  }
  //
  // Muon Trigger
  //
  if(iPWG1perfMuonTrig) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskMTRchamberEfficiency.C");
    AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
    if(!taskMuonTrig) { 
      Error("RunPerformanceTrain","AliAnalysisTaskTrigChEff not created!");
      return;
    }
    mgr->AddTask(taskMuonTrig);
  } 
  else {
    Printf("RunPerformanceTrain: AliAnalysisTaskTrigChEff - EXCLUDED!");
  }
  //
  // Muon Efficiency
  //
  if(iPWG1perfMuonEff) {
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
  AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AliAnalysisTaskMuonTrackingEff();
  if(!taskMuonTrackEff) { 
     Error("RunPerformanceTrain","AliAnalysisTaskMuonTrackingEff not created!");
     return;
  }
  mgr->AddTask(taskMuonTrackEff);
  } 
  else {
    Printf("RunPerformanceTrain: Muon Efficiency - EXCLUDED!");
  }
  //
  // TOF performance
  //
  if(iPWG1perfTOF) {
  //
  } 
  else {
    Printf("RunPerformanceTrain: TOF - EXCLUDED!");
  }
  //
  // PWG1 Primary Vertex
  //
  if(iPWG1perfPrimVertex) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD *taskPrimVertex = AddTaskVertexESD();
    if(!taskPrimVertex) { 
      Error("RunPerformanceTrain","AliAnalysisTaskVertexESD not created!");
      return;
    }
  } 
  else {
    Printf("RunPerformanceTrain: AliAnalysisTaskVertexESD - EXCLUDED!");
  }
  //
  // PWG1 V0 QA
  //
  if (iPWG1v0QA) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskV0QA.C");
    AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(bUseMCInfo);
    if(!taskv0QA) {
      Error("RunPerformanceTrain","AliAnalysisTaskV0QA not created!");
      return;
    }    
  }
  else {
    Printf("RunPerformanceTrain: AliAnalysisTaskV0QA - EXCLUDED!");
  }

  //
  // Disable debug printouts
  //
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  //mgr->StartAnalysis("local",chain, nEvents, firstEvent);
  mgr->StartAnalysis("local",chain);
}

