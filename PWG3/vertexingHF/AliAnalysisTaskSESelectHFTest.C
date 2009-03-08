void AliAnalysisTaskSESelectHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour selection
  // and creation of a stand-alone AOD
  // A.Dainese, andrea.dainese@lnl.infn.it
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
  //inputHandler->AddFriend("AliAOD.VertexingHF.root");
  mgr->SetInputEventHandler(inputHandler);

  // Output 
  AliAODHandler *aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.VertexingHF.sa.root");
  aodHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodHandler);

  
  // Aanalysis task    
  AliAnalysisTaskSESelectHF *hfTask = new AliAnalysisTaskSESelectHF("SelectHFAnalysis");
  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);
  
  //
  // Create containers for input/output
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());
  /*
  // before v4-17-Release
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							    AliAnalysisManager::kOutputContainer, 
							    "default");
  mgr->ConnectInput(hfTask,0,cinput1);
  mgr->ConnectOutput(hfTask,0,coutput1);
  */

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
