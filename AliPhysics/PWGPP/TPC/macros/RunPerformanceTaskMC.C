// Macro to run TPC MC performance
//
// The basic tracking functionality is tested 
// by using MC track references (exact points). 
// 
// Track references are propagate in B-field and
// material by using:
//   
// AliTracker::PropagateTrackToBxByBz()
// AliExternalTrackParam::PropagateToBxByBz()
//
// 
// Available test components:
//
// AliPerformanceRes (propagation TPCin(ref) -> DCA(particle) by using AliTracker::PropagateTrackToBxByBz() and comparison at DCA)
// AliPerformanceResTCPInner (propagation TPCout(ref) -> TPCin(ref) and comparison at TPCin)
// AliPerformanceResTPCOuter (propagation TPCin(ref) -> TPCout(ref) by using AliExternalTrackParam::PropagateToBxByBz() and comparison at TPCout) 
//

/*
  //1. Run locally e.g.

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/LoadMyLibs.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list_flatP_JB.txt", 500, 0);
  chain->Lookup();

  // Geometry (need for the track propagation through material)
  //AliGeomManager::LoadGeometry("/lustre/alice/jacek/sim/HEADJB/flatPt_uniB/0/geometry.root");

  // set magnetic field
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/RunPerformanceTaskMC.C");
  RunPerformanceTaskMC(chain, kTRUE, kFALSE, kFALSE, 0);

  //2. Run on PROOF Lite e.g.

  TProof::Open(""); 

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/d/alice11/jacek/alice/x86_64/AliRoot/trunkJB");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list_flatP_JB.txt", 400, 0);
  chain->Lookup();

  // set magnetic field
  // the best is to create macro MagField.C with the line: 
  // TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
  gProof->Exec("gROOT->Macro(\"MagField.C\")");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/RunPerformanceTaskMC.C");
  RunPerformanceTaskMC(chain, kTRUE, kTRUE, kTRUE,0);

  //3. Run only on static PROOF at GSI e.g.

  TProof::Reset("jacek@lxgrid5.gsi.de");
  TProofMgr * proofmgr = TProof::Mgr("jacek@lxgrid5.gsi.de");
  //proofmgr->SetROOTVersion("523-04");
  TProof * proof = proofmgr->CreateSession();
  proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEADJB/");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("flat_JB.txt", 50, 0);
  chain->Lookup();

  // Geometry (need for the track propagation through material)
  //AliGeomManager::LoadGeometry("/lustre/alice/local/TRDdata/SIM/P-Flat/TRUNK/test/RUN0/geometry.root");

  // set magnetic field
  gProof->Exec("gROOT->ProcessLine(\"TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\", 1., 1., AliMagF::k5kG))\")");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/RunPerformanceTaskMC.C");
  RunPerformanceTaskMC(chain, kTRUE, kTRUE, kTRUE);

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
void RunPerformanceTaskMC(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE,  Bool_t bProof=kTRUE, Int_t debugStreamLevel=0)
{
  if(!chain) 
  {
    AliDebug(AliLog::kError, "ERROR: No input chain available");
    return;
  }

  //
  // Create global cuts objects 
  //

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts("pRecInfoCuts"); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetPtRange(0.15,1.e10);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts("pMCInfoCuts");
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create performance objects and set cuts 
  //
  const Int_t kTPC = 0; const Int_t kTPCITS = 1; const Int_t kConstrained = 2; const Int_t kTPCInner = 3; const Int_t kTPCOuter = 4;

  //
  // MC 
  //
  AliPerformanceMC *pCompMC0 = new AliPerformanceMC("AliPerformanceMC","AliPerformanceMC",kTPC,kFALSE); 
  if(!pCompMC0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceMC object");
  }
  pCompMC0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompMC0->SetAliMCInfoCuts(pMCInfoCuts);
  //pCompMC0->SetStreamLevel(debugStreamLevel);

  AliPerformanceMC *pCompMC3 = new AliPerformanceMC("AliPerformanceMCTPCInner","AliPerformanceMCTPCInner",kTPCInner,kFALSE); 
  if(!pCompMC3) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceMCInnerTPC object");
  }
  pCompMC3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompMC3->SetAliMCInfoCuts(pMCInfoCuts);
  //pCompMC3->SetStreamLevel(debugStreamLevel);

  AliPerformanceMC *pCompMC4 = new AliPerformanceMC("AliPerformanceMCTPCOuter","AliPerformanceMCTPCOuter",kTPCOuter,kFALSE); 
  if(!pCompMC4) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceMCTPCOuter object");
  }
  pCompMC4->SetAliRecInfoCuts(pRecInfoCuts);
  pCompMC4->SetAliMCInfoCuts(pMCInfoCuts);
  //pCompMC4->SetStreamLevel(debugStreamLevel);


  // Swtich off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager;

  // Create task
  AliPerformanceTask *task = new AliPerformanceTask("PerformanceMC","TPC Performance");
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);
  task->AddPerformanceObject( pCompMC0 );
  task->AddPerformanceObject( pCompMC3 );
  task->AddPerformanceObject( pCompMC4 );

  // Add task
  mgr->AddTask(task);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  if(bUseESDfriend) esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);

  if(bUseMCInfo) {
  // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(handler);
  }

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

