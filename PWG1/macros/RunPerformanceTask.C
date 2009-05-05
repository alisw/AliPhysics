// Macro to run TPC performance

/*
  //1. Run locally e.g.

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/LoadMyLibs.C");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esd_v4-16-Rev-08-grid.txt", 10, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kFALSE);

  //2. Run on PROOF Lite e.g.

  TProof::Open(""); 

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esd_v4-16-Rev-08-grid.txt", 500, 0);
  //TChain* chain = CreateESDChain("esd_TRUNK_flat_ideal_geom.txt", 50, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kTRUE); 


  //3. Run only on static PROOF at GSI e.g.

  TProof::Reset("jacek@lxgrid5.gsi.de");
  TProofMgr * proofmgr = TProof::Mgr("jacek@lxgrid5.gsi.de");
  proofmgr->SetROOTVersion("523-02");
  TProof * proof = proofmgr->CreateSession();
  proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esd_v4-16-Rev-08-grid.txt", 200, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kTRUE); 

  //4. Make final spectra and store them in the
  // output folder and generate control pictures e.g.

  TFile f("TPC.Performance.root");
  AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceResTPCInner");
  compObjRes->Analyse();
  compObjRes->GetAnalysisFolder()->ls("*");
  compObjRes->PrintHisto(kTRUE,"PerformanceResTPCInnerQA.ps");
  TFile fout("AnalysedResTPCInner.root","recreate");
  compObjRes->GetAnalysisFolder()->Write();
  fout.Close();
  f.Close();


*/

//_____________________________________________________________________________
void RunPerformanceTask(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bProof=kTRUE)
{
  if(!chain) 
  {
    AliDebug(AliLog::kError, "ERROR: No input chain available");
    return;
  }

  // set Proof
  if(bProof) { 

    /*
    // ONLY FOR GSI
    TProof::Reset("jacek@lxgrid5.gsi.de");
    TProofMgr * proofmgr = TProof::Mgr("jacek@lxgrid5.gsi.de");
    proofmgr->SetROOTVersion("523-02");
    TProof * proof = proofmgr->CreateSession();
    proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);
    */

    /*
    //cout << "*** START PROOF Lite SESSION ***" << endl;
    TProof::Open(""); 
    gROOT->LoadMacro("ProofEnableAliRoot.C");
    ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");
    */

  }
  // set magnetic field
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  //
  // Create global cuts objects 
  //

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create performance objects and set cuts 
  //
  const Int_t kTPC = 0; const Int_t kTPCITS = 1; const Int_t kConstrained = 2; const Int_t kTPCInner = 3;

  //
  // Resolution
  //
  AliPerformanceRes *pCompRes0 = new AliPerformanceRes("AliPerformanceRes","AliPerformanceRes",kTPC,kFALSE); 
  if(!pCompRes0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceRes object");
  }
  pCompRes0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);

  AliPerformanceRes *pCompRes3 = new AliPerformanceRes("AliPerformanceResTPCInner","AliPerformanceResTPCInner",kTPCInner,kFALSE); 
  if(!pCompRes3) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceResInnerTPC object");
  }
  pCompRes3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes3->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // Efficiency
  //
  AliPerformanceEff *pCompEff0 = new AliPerformanceEff("AliPerformanceEff","AliPerformanceEff",kTPC,kFALSE); 
  if(!pCompEff0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceEff object");
  }
  pCompEff0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // dEdx
  //
  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceDEdxTPCInner object");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // DCA
  //
  AliPerformanceDCA *pCompDCA0 = new AliPerformanceDCA("AliPerformanceDCA","AliPerformanceDCA",kTPC,kFALSE); 
  if(!pCompDCA0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceDCA object");
  }
  pCompDCA0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC performance
  //
  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE); 
  if(!pCompTPC0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceTPC object");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
 
  // Swtich off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager;

  // Create task
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");
  if (bUseMCInfo) task->SetUseMCInfo(kTRUE);
  task->AddPerformanceObject( pCompRes0 );
  task->AddPerformanceObject( pCompRes3 );
  task->AddPerformanceObject( pCompEff0 );
  task->AddPerformanceObject( pCompDEdx3 );
  task->AddPerformanceObject( pCompDCA0 );
  task->AddPerformanceObject( pCompTPC0 );

  // Add task
  mgr->AddTask(task);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  //esdH->SetInactiveBranches("*");
  mgr->SetInputEventHandler(esdH);

  if(bUseMCInfo) {
  // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    handler->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Create input chain
  //gROOT->LoadMacro("CreateESDChain.C");
  //TChain* chain = CreateESDChain(fileList, NumberOfFiles, fromFile);
  //if(!chain) {
  //  printf("ERROR: chain cannot be created\n");
  //  return;
  //}

  // Create containers for input
  //AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  // Create containers for output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if(bProof) mgr->StartAnalysis("proof",chain);
  else mgr->StartAnalysis("local",chain);
}

