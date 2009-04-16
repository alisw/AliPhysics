void AliAnalysisTaskSEBtoJPSItoEleTest() 
{
  //
  // Test macro for the AliAnalysisTaskSEBtoJPSItoEle 
  // starting from AliAOD.root file with HF + Like Sign candidates.
  // C.Di Giglio, carmelo.digiglio@ba.infn.it
  //

  TString mode="local"; // otherwise, "grid" 
  //TString mode="grid"; // needs definition of a proper Grid handler  

  // load proper libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3vertexingHF.so");
 
  TChain *chain = 0;

  if(mode=="local") {
    // Local files 
    TString treeName,fileName;
    treeName="aodTree";
    fileName="AliAOD.root";
    chain = new TChain(treeName.Data());
    chain->Add(fileName.Data());
  } else if (mode=="grid") {
    // Fetch files with AliEn :
    const char *collectionfile = "Collection.xml";
    TGrid::Connect("alien://") ;
    TAlienCollection *coll   = TAlienCollection::Open(collectionfile);
    chain = new TChain("aodTree");
    while(coll->Next()) chain->Add(coll->GetTURL(""));
  } else {
    printf("ERROR: mode has to be \"local\" or \"grid\" \n");
    return;
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler();
  inputHandler->AddFriend("AliAOD.VertexingHF.root");
  mgr->SetInputEventHandler(inputHandler);

  // Cdf unbinned log-likelihood fit analysis task    
  AliAnalysisTaskSEBtoJPSItoEle *hfTask = new AliAnalysisTaskSEBtoJPSItoEle("CdfFitAnalysis");
  hfTask->SetDebugLevel(2);

  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  // before v4-17-Release
  /*AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(),
                                                           AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            "default");

  mgr->ConnectInput(hfTask,0,cinput1);
  mgr->ConnectOutput(hfTask,0,coutput1);*/

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "CdfFit.root");

  mgr->ConnectOutput(hfTask,1,coutput);

  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if(!mgr->InitAnalysis()) return;

  mgr->PrintStatus();

  mgr->StartAnalysis(mode.Data(),chain);

  return;

}
