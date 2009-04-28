void RunAnalysisAODVertexingHF()
{
  //
  // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
  // It has the structure of a Analysis Train:
  // - in this macro, change things related to running mode
  //   and input preparation 
  // - add your task using a AddTaskXXX macro 
  //
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //

  TString analysisMode = "local"; // "local" or "grid"
  TString inputMode    = "list"; // "list" or "xml"


  Bool_t useParFiles=kFALSE;


  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/MakeAODInputChain.C");
  LoadLibraries(useParFiles);

  if(analysisMode=="grid") TGrid::Connect("alien:",0,0,"t");




  //-------------------------------------------------------------------
  // Prepare input chain
  TChain *chainAOD = 0;
  
  if(inputMode=="list") {
    // Local files
    chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
    //chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/analysis/out_lhcw/290001/",2,2);
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
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks    
  //
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskCompareHF.C");
  AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskSelectHF.C");
  //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();
  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskBkgLikeSign.C");
  //AliAnalysisTaskSEBkgLikeSignJPSI *lsTask = AddTaskBkgLikeSign();
  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskBtoJPSItoEle.C");
  //AliAnalysisTaskSEBtoJPSItoEle *jpsiTask = AddTaskBtoJPSItoEle();
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskCF.C");
  AliCFHeavyFlavourTask *cfTask = AddTaskCF();
  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskCFMultiVar.C");
  //AliCFHeavyFlavourTaskMultiVar *cfmvTask = AddTaskCFMultiVar();


  //-------------------------------------------------------------------

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
