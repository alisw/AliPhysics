void AliAnalysisTaskMEVertexingHFTest()
{
  //
  // Test macro for the AliAnalysisTaskME for heavy-flavour event mixing
  // r.romita@gsi.de
  //

  Bool_t useParFiles=kFALSE;
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  // Local files 
  

  TChain* chain = new TChain("aodTree");
  Char_t fileName[100];
  sprintf(fileName,"AliAODs.root");
  chain->Add(fileName);
  
  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  
  // Input Handler
  AliMultiEventInputHandler *inputHandler = new AliMultiEventInputHandler(4,1);
  AliEventPoolOTF* pool = new AliEventPoolOTF("event pool", "AOD");
  // apply selections
  pool->SetMultiplicityBin(0, 100, 2);
  pool->SetZVertexBinning(-20., 20., 2);
  pool->Init();
  //set tag directory
  Char_t tagDir[100];
  sprintf(tagDir,".");
  pool->SetTagDirectory(tagDir);
  mgr->SetInputEventHandler(inputHandler);
  mgr->SetEventPool(pool);
  inputHandler->SetEventPool(pool);
  
  // Output 
  AliAODHandler *aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.VertexingHF.root");
  aodHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodHandler);
  
  gROOT->LoadMacro("AddTaskMixing.C");
  AliAnalysisTaskMEVertexingHF *hfTask = AddTaskHFMixing();
  
  
  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if(!mgr->InitAnalysis()) return;

  mgr->PrintStatus();

  TStopwatch watch;
  watch.Start();
  mgr->StartAnalysis("mix",chain, 1000);
  watch.Stop();
  watch.Print();
  delete mgr;

  return;
}
