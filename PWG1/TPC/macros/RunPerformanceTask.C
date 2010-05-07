// Macro to run TPC performance task (locally, proof).
//
// By default 8 performance components are added to 
// the task. Look inside AddTaskPerformanceTPC.C
//

/*
 
  //1. Run locally e.g.

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/LoadMyLibs.C");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esds_test.txt",10, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kFALSE, kTRUE, kFALSE);

  //2. Run on PROOF Lite e.g.

  TProof::Open(""); 

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/trunk");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list_flatP_JB.txt",20, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kTRUE, kTRUE);


  //3. Run only on PoD at GSI e.g.
  TProofMgr * proofmgr = TProof::Mgr("lxialpod2.gsi.de:21001");
  TProof * proof = proofmgr->CreateSession();
  proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/d/alice11/jacek/alice/x86_64/AliRoot/trunk");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("../input/ffprod_v4-17-Rev-19_900kPythia6D6T.list", 200, 0);
  chain->Lookup();

  //gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/RunPerformanceTask.C");
  gROOT->LoadMacro("/d/alice11/jacek/alice/TPC/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kFALSE, kTRUE, kTRUE); 

  //4. Make final spectra and store them in the
  // output folder and generate control pictures e.g.

  TFile f("TPC.Performance.root");
  AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceResTPCOuter");
  compObjRes->Analyse();
  compObjRes->GetAnalysisFolder()->ls("*");
  // create pictures
  compObjRes->PrintHisto(kTRUE,"PerformanceResTPCOuterQA.ps");
  // store output QA histograms in file
  TFile fout("PerformanceResTPCOuterQAHisto.root","recreate");
  compObjRes->GetAnalysisFolder()->Write();
  fout.Close();
  f.Close();

*/

//_____________________________________________________________________________
void RunPerformanceTask(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE,  Bool_t bProof=kTRUE)
{
  if(!chain) 
  {
    AliDebug(AliLog::kError, "ERROR: No input chain available");
    return;
  }
  //
  // Swtich off all AliInfo (too much output!)
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //
  // Create analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager;
  if(!mgr) { 
    Error("runTPCQA","AliAnalysisManager not set!");
    return;
  }

  //
  // Set ESD input handler
  //
  AliESDInputHandler* esdH = new AliESDInputHandler;
  if(!esdH) { 
    Error("runTPCQA","AliESDInputHandler not created!");
    return;
  }
  if(bUseESDfriend) esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);

  //
  // Set MC input handler
  //
  if(bUseMCInfo) {
    AliMCEventHandler* mcH = new AliMCEventHandler;
    if(!esdH) { 
      Error("runTPCQA","AliMCEventHandler not created!");
      return;
    }
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH);
  }
  //
  // Add task to AliAnalysisManager
  //
  //gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPC.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/AddTaskPerformanceTPCQA.C");
  AliPerformanceTask *tpcQA = AddTaskPerformanceTPCQA(bUseMCInfo,bUseESDfriend);
  if(!tpcQA) { 
      Error("runTPCQA","TaskPerformanceTPC not created!");
      return;
  }

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if(bProof) mgr->StartAnalysis("proof",chain);
  else mgr->StartAnalysis("local",chain);
}

