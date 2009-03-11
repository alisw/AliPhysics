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
