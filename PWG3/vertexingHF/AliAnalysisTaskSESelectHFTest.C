void AliAnalysisTaskSESelectHFTest()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour selection
  // and creation of a stand-alone AOD
  // A.Dainese, andrea.dainese@lnl.infn.it
  //

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
  Double_t cutsD0[9]=
    // cutsD0[0] = inv. mass half width [GeV]
    // cutsD0[1] = dca [cm]
    // cutsD0[2] = cosThetaStar
    // cutsD0[3] = pTK [GeV/c]
    // cutsD0[4] = pTPi [GeV/c]
    // cutsD0[5] = d0K [cm]   upper limit!
    // cutsD0[6] = d0Pi [cm]  upper limit!
    // cutsD0[7] = d0d0 [cm^2]
    // cutsD0[8] = cosThetaPoint
                     {1,
		      100000.,
		      1.1,
		      0.,
		      0.,
		      100000.,
		      100000.,
		      100000000.,
		      -0.9}; 
  hfTask->SetD0toKpiCuts(cutsD0);
  
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
