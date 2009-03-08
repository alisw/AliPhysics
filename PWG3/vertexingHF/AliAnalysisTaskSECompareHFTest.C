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
  LoadLibraries(useParFiles);

  if(analysisMode=="grid") TGrid::Connect("alien:",0,0,"t");




  TChain *chainAOD = 0;
  TChain *chainAODfriend = 0;
  
  if(inputMode=="list") {
    // Local files
    chainAOD = new TChain("aodTree");
    chainAODfriend = new TChain("aodTree");
    // set the path to the files (can be local or on alien)
    chainAOD->Add(      "alien:///alice/cern.ch/user/r/rbala/analysis/out_lhcw/290001/2/AliAOD.root");
    chainAODfriend->Add("alien:///alice/cern.ch/user/r/rbala/analysis/out_lhcw/290001/2/AliAOD.VertexingHF.root");
    // ... add more if needed
  } else if(inputMode=="xml") {
    // xml: need to check the 1-to-1 correspondence
    TString collectionfileAOD       = "collection_aod.xml";
    TString collectionfileAODfriend = "collection_aodHF.xml";
    TAlienCollection *collectionAOD       = TAlienCollection::Open(collectionfileAOD.Data());
    TAlienCollection *collectionAODfriend = TAlienCollection::Open(collectionfileAODfriend.Data());
    chainAOD = new TChain("aodTree");
    chainAODfriend = new TChain("aodTree");
    while(collectionAOD->Next())       chainAOD->Add(collectionAOD->GetTURL(""));
    while(collectionAODfriend->Next()) chainAODfriend->Add(collectionAODfriend->GetTURL(""));
  }

  // attach the friend chain
  chainAOD->AddFriend(chainAODfriend);

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
