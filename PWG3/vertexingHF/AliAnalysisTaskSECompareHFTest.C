void AliAnalysisTaskSECompareHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour candidates
  // association with MC truth (using MC info in AOD)
  // A.Dainese, andrea.dainese@lnl.infn.it
  //
  TString mode="local"; // otherwise, "grid"

  Bool_t useParFiles=kFALSE;

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  TChain *chainAOD = 0;
  TChain *chainAODfriend = 0;

  chainAOD = new TChain("aodTree");
  chainAODfriend = new TChain("aodTree");

  if(mode=="local") {
    // Local files 
    chainAOD->Add("./AliAOD.root");
    chainAODfriend->Add("./AliAOD.VertexingHF.root");
    // ... add more if needed
  } else { // grid (NOT TESTED YET)
    //Fetch files with AliEn
    // xml: need to check the 1-to-1 correspondence
    const char *collectionfileAOD = "CollectionAOD.xml";
    const char *collectionfileAODfriend = "CollectionAODfriend.xml";
    TGrid::Connect("alien://");
    //Create an AliTagAnalysis Object and chain the tags
    AliTagAnalysis   *tagAna = new AliTagAnalysis();
    chainAOD = tagAna->GetChainFromCollection(collectionfileAOD,"aodTree");
    chainAODfriend = tagAna->GetChainFromCollection(collectionfileAODfriend,"aodTree");
  }

  // attach the friend chain
  chainAOD->AddFriend(chainAODfriend);

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  //inputHandler->AddFriend("AliAOD.VertexingHF.root");//should not be needed
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
  mgr->ConnectInput(hfTask,0,cinput);
  mgr->ConnectOutput(hfTask,1,coutput);

  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  if(mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis(mode.Data(),chainAOD);
  }

  return;
}
