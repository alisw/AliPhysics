// Macro to run TPC performance task (locally, proof).
//
// By default 8 performance components are added to 
// the task. Look inside AddTaskPerformanceTPC.C
//

/*
 
  //1. Run locally e.g.

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/LoadMyLibs.C");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esds_test.txt",10, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kFALSE, kTRUE, kFALSE);

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
void RunPerformanceTask(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE,  Bool_t /*legacyMode*/=kFALSE)
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
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskPerformanceTPC.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskPerformanceTPCQA.C");
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

  mgr->StartAnalysis("local",chain);
}

