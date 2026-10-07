void Load(const char* taskName, Bool_t debug)
{
  TString compileTaskName;
  compileTaskName.Form("%s.cxx++", taskName);
  if (debug)
    compileTaskName += "g";

  gROOT->Macro(compileTaskName);

  // Enable debug printouts
  if (debug)
  {
    AliLog::SetClassDebugLevel(taskName, AliLog::kDebug+2);
  }
  else
    AliLog::SetClassDebugLevel(taskName, AliLog::kWarning);
}

void run(const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t inputMode = kFALSE, const char* option = "")
{
  // inputMode option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset
  //
  // option is passed to the task(s)

  if (nRuns < 0)
    nRuns = 1234567890;

  {
    gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/"); 
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
    gSystem->Load("libPWG0dep");
  }
  
  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades AliESDTZERO MuonTracks TrdTracks CaloClusters");
  mgr->SetInputEventHandler(esdH);

  cInput = mgr->GetCommonInputContainer();
  
  Load("AliEventStatsTask", aDebug);
  TString optStr(option);
  
  // remove SAVE option if set
  Bool_t save = kFALSE;
  if (optStr.Contains("SAVE"))
  {
    optStr = optStr(0,optStr.Index("SAVE")) + optStr(optStr.Index("SAVE")+4, optStr.Length());
    save = kTRUE;
  }
  
  task = new AliEventStatsTask(optStr);
  physicsSelection = new AliPhysicsSelection;
  if (aDebug)
    AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);
  task->SetPhysicsSelection(physicsSelection);
  //AliBackgroundSelection* background = new AliBackgroundSelection("AliBackgroundSelection", "AliBackgroundSelection");
  //physicsSelection->AddBackgroundIdentification(background);
  
  mgr->AddTask(task);

  // Attach input
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 1, cOutput);

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();

  if (inputMode == -1)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    TFile::Open("dataset.root");
    ds = (TFileCollection*) gFile->Get("dataset");
    chain = CreateChainFromDataSet(ds);
    mgr->StartAnalysis("local", chain, nRuns, offset);
  }
  else
  {
    // Create chain of input files
    gROOT->LoadMacro("../CreateESDChain.C");

    chain = CreateESDChain(data, nRuns, offset);
    //chain = CreateChain("TE", data, nRuns, offset);

    mgr->StartAnalysis("local", chain);
  }

}
