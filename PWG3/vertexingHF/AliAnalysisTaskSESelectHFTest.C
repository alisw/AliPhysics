void AliAnalysisTaskSESelectHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour selection
  // and creation of a stand-alone AOD
  // A.Dainese, andrea.dainese@lnl.infn.it
  //

  Bool_t useParFiles=kFALSE;

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);


  // Local files 
  TChain *chain = new TChain("aodTree");
  chain->Add("./AliAOD.root");

  // or:
  /*
  //Fetch files with AliEn :
  const char *collectionfile = "CollectionTags.xml";
  TGrid::Connect("alien://") ;
  //Create an AliRunTagCuts and an AliEventTagCuts Object and impose some selection criteria
  AliRunTagCuts      *runCuts   = new AliRunTagCuts();
  AliEventTagCuts    *eventCuts = new AliEventTagCuts();
  AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts();
  eventCuts->SetMultiplicityRange(0,20000);
  //Create an AliTagAnalysis Object and chain the tags
  AliTagAnalysis   *tagAna = new AliTagAnalysis();
  tagAna->SetType("AOD");
  TAlienCollection *coll   = TAlienCollection::Open(collectionfile);
  TGridResult      *tagResult = coll->GetGridResult("",0,0);
  tagResult->Print();
  tagAna->ChainGridTags(tagResult);
  //Create a new aod chain and assign the chain that is returned by querying the tags
  TChain* chain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts);
  */


  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  inputHandler->AddFriend("AliAOD.VertexingHF.root");
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
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if(mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
    //mgr->StartAnalysis("grid",chain);
  }

  return;
}
