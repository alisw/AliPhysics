void AliAnalysisTaskSECompareHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour candidates
  // association with MC truth (using MC info in AOD)
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //

  TString analysisMode = "grid"; // "local" or "grid"
  TString inputMode    = "list"; // "list" or "xml"


  Bool_t useParFiles=kFALSE;


  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/MakeAODInputChain.C");
  LoadLibraries(useParFiles);

  if(analysisMode=="grid") TGrid::Connect("alien:",0,0,"t");




  TChain *chainAOD = 0;
  
  if(inputMode=="list") {
    // Local files
    //chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
    chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/analysis/out_lhcw/290001/",2,2);
  } else if(inputMode=="xml") {
    // xml
    chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputHandler);

  
  // Aanalysis task    
  AliAnalysisTaskSECompareHF *hfTask = new AliAnalysisTaskSECompareHF("CompareHFAnalysis");
  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cinput",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "CmpHF.root");
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(hfTask,1,coutput);

  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  if(mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chainAOD);
  }

  return;
}
