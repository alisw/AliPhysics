void AliAnalysisTaskSEVertexingHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour vertexing
  // A.Dainese, andrea.dainese@lnl.infn.it
  //

  Bool_t inputAOD=kFALSE; // otherwise, ESD
  Bool_t createAOD=kTRUE; // kTRUE: create AOD and use it as input to vertexing
                          // kFALSE: use ESD as input to vertexing
  Bool_t writeKineToAOD = kFALSE;
  TString mode="local"; // otherwise, "grid" 
  Bool_t useParFiles=kFALSE;
  Bool_t doCentrality=kTRUE;

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");
  LoadLibraries(useParFiles);
  gSystem->Load("libPWG3muon");
  TChain *chain = 0;

  if(mode=="local") {
    // Local files 
    TString treeName,fileName;
    if(inputAOD) {
      treeName="aodTree"; 
      fileName="AliAOD.root";
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
  const char* deltaAODfname="AliAOD.VertexingHF.root";
  if(createAOD) {
    aodHandler->SetOutputFileName("AliAOD.root");
  } else {
    aodHandler->SetOutputFileName(deltaAODfname);
    aodHandler->SetAODExtensionMode();
  }
  mgr->SetOutputEventHandler(aodHandler);
  mgr->RegisterExtraFile(deltaAODfname);  

  if(!inputAOD && createAOD) {
    // MC Truth
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    if(writeKineToAOD) mgr->SetMCtruthEventHandler(mcHandler);
    AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
    if(writeKineToAOD) mgr->AddTask(kinefilter);  
    // Centrality
    if(doCentrality){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    }

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
    AliAnalysisTaskESDfilter *filter = AddTaskESDFilter(writeKineToAOD);
   
  }

  // Vertexing analysis task    
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/AddTaskVertexingHF.C");
  AliAnalysisTaskSEVertexingHF *hfTask = AddTaskVertexingHF(deltaAODfname);
  
  
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
