void AliAnalysisTaskSEVertexingHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour vertexing
  // A.Dainese, andrea.dainese@lnl.infn.it
  //
  Bool_t inputAOD=kTRUE; // otherwise, ESD


  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3vertexingHF.so");


  // Local files 
  TString treeName,fileName;
  if(inputAOD) {
    treeName="aodTree"; 
    fileName="AliAOD.root";
  } else {
    treeName="esdTree"; 
    fileName="AliESDs.root";
  }
  TChain *chain = new TChain(treeName.Data());
  chain->Add(fileName.Data());

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
  tagAna->SetType("ESD");
  TAlienCollection *coll   = TAlienCollection::Open(collectionfile);
  TGridResult      *tagResult = coll->GetGridResult("",0,0);
  tagResult->Print();
  tagAna->ChainGridTags(tagResult);
  //Create a new esd chain and assign the chain that is returned by querying the tags
  TChain* chain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts);
  */


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
  AliAnalysisTaskSEVertexingHF *hfTask = new AliAnalysisTaskSEVertexingHF("VertexingHFAnalysis");
  hfTask->SetDebugLevel(2);
  
  mgr->AddTask(hfTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							    AliAnalysisManager::kOutputContainer, 
							    "default");
  mgr->ConnectInput(hfTask,0,cinput1);
  mgr->ConnectOutput(hfTask,0,coutput1);

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
