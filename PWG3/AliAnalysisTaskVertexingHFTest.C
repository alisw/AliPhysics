void AliAnalysisTaskVertexingHFTest() {

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISRL");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");
  //This file can cause problems :
  if (!gSystem->AccessPathName("$ALICE_ROOT/ANALYSIS/AliAnalysisSelector_cxx.so",kFileExists)) {
    printf("File $ALICE_ROOT/ANALYSIS/AliAnalysisSelector_cxx.so exists and can cause problems, delete it if you can...");
    return;
  }

 
  //Run over local files :
  TChain *chain= new TChain("esdTree");
  chain->Add("AliESDs.root"); // put path to your files here
  
  // OR :
 
  /*
  //Fetch files with AliEn :
  const char *collectionfile = "essai6000CollectionTags1.xml";
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
 
  //Temporary solution to avoid memory leaks : 
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*CaloClusters*",0);

  //Create tasks
  AliAnalysisManager *analManager = new AliAnalysisManager("myAnalysisManager");
  //analManager->SetDebugLevel(10);
  AliAnalysisTaskVertexingHF *task1 = new AliAnalysisTaskVertexingHF("myTask");
  analManager->AddTask(task1);
  //Create containers for input/output
  AliAnalysisDataContainer *cinput1 = analManager->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = analManager->CreateContainer("tree1", TTree::Class(),AliAnalysisManager::kOutputContainer,"HFtrees.root");
  AliAnalysisDataContainer *coutput2 = analManager->CreateContainer("tree2", TTree::Class(),AliAnalysisManager::kOutputContainer,"HFtrees.root");
  AliAnalysisDataContainer *coutput3 = analManager->CreateContainer("tree3", TTree::Class(),AliAnalysisManager::kOutputContainer,"HFtrees.root");
  AliAnalysisDataContainer *coutput4 = analManager->CreateContainer("tree4", TTree::Class(),AliAnalysisManager::kOutputContainer,"HFtrees.root");
  analManager->ConnectInput(task1,0,cinput1);
  analManager->ConnectOutput(task1,0,coutput1);
  analManager->ConnectOutput(task1,1,coutput2);
  analManager->ConnectOutput(task1,2,coutput3);
  analManager->ConnectOutput(task1,3,coutput4);
  cinput1->SetData(chain);
 
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if (analManager->InitAnalysis()) {
    analManager->PrintStatus();
    //analManager->StartAnalysis("grid",chain);
    analManager->StartAnalysis("local",chain);
  }
 
  return;
}



