void AliAnalysisTaskSEVertexingHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour vertexing
  // A.Dainese, andrea.dainese@lnl.infn.it
  //

  Bool_t inputAOD=kTRUE; // otherwise, ESD
  TString mode="local"; // otherwise, "grid" 
  Bool_t useParFiles=kFALSE;

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  TChain *chain = 0;

  if(mode=="local") {
    // Local files 
    TString treeName,fileName;
    if(inputAOD) {
      treeName="aodTree"; 
      fileName="AliAODs.root";
    } else {
      treeName="esdTree"; 
      fileName="AliESDs.root";
    }
    chain = new TChain(treeName.Data());
    chain->Add(fileName.Data());

  } else if (mode=="grid") {
    //Fetch files with AliEn :
    const char *collectionfile = "Collection.xml";
    TGrid::Connect("alien://") ;
    TAlienCollection *coll   = TAlienCollection::Open(collectionfile);
    if(inputAOD) { // input AOD
      chain = new TChain("aodTree");
      while(coll->Next()) chain->Add(coll->GetTURL(""));
    } else { // input ESD
      //Create an AliRunTagCuts and an AliEventTagCuts Object and impose some selection criteria
      AliRunTagCuts      *runCuts   = new AliRunTagCuts();
      AliEventTagCuts    *eventCuts = new AliEventTagCuts();
      AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts();
      AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts();
      eventCuts->SetMultiplicityRange(0,20000);
      //Create an AliTagAnalysis Object and chain the tags
      AliTagAnalysis   *tagAna = new AliTagAnalysis();
      tagAna->SetType("ESD");
      TGridResult      *tagResult = coll->GetGridResult("",0,0);
      tagResult->Print();
      tagAna->ChainGridTags(tagResult);
      //Create a new esd chain and assign the chain that is returned by querying the tags
      chain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts);
    }
  } else {
    printf("ERROR: mode has to be \"local\" or \"grid\" \n");
    return;
  }


  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  
  // Input Handler
  AliInputEventHandler *inputHandler = 0;
  if(inputAOD) {
    inputHandler = new AliAODInputHandler();
  } else {
    inputHandler = new AliESDInputHandler();
  }
  mgr->SetInputEventHandler(inputHandler);
  
  // Output 
  AliAODHandler *aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.VertexingHF.root");
  aodHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodHandler);
  
  // Vertexing analysis task    
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskVertexingHF.C");
  AliAnalysisTaskSEVertexingHF *hfTask = AddTaskVertexingHF();
  
  
  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if(!mgr->InitAnalysis()) return;

  mgr->PrintStatus();

  TStopwatch watch;
  watch.Start();
  mgr->StartAnalysis(mode.Data(),chain);
  watch.Stop();
  watch.Print();

  return;
}
