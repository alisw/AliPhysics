//
// Macro to run performance QA train
// locally. The TPC performance task is attached.
// 
//
//13.10.2009 -  J.Otwinowski@gsi.de

//_____________________________________________________________________________
void RunPerformanceTrain(Char_t *list="esdList.txt", Int_t nFiles=20, Int_t fistFile=0, Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE)
{
  //
  // Swtich off all AliInfo (too much output!)
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //
  // Create input ESD chain
  //
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(list,nFiles,fistFile);
  if(!chain) { 
    Error("RunPerformanceTrain","ESD chain not created!");
    return;
  }
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

  //
  // Set MC input handler
  //
  if(bUseMCInfo) {
    AliMCEventHandler* mcH = new AliMCEventHandler;
    if(!esdH) { 
      Error("RunPerformanceTrain","AliMCEventHandler not created!");
      return;
    }
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH);
  }
  //
  // Add task to AliAnalysisManager
  //
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPC.C");
  AliPerformanceTask *tpcQA = AddTaskPerformanceTPC(bUseMCInfo,bUseESDfriend);
  if(!tpcQA) { 
      Error("RunPerformanceTrain","TaskPerformanceTPC not created!");
      return;
  }

  //
  // Disable debug printouts
  //
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

