void RunAlignmentDataFilterITS() {
  //
  // Macro to extract AliTrackPoints for ITS
  // A.Dainese, andrea.dainese@pd.infn.it
  //

  Int_t nentries=1234567890;
  Int_t firstentry=0;

  // Load geometry file
  AliGeomManager::LoadGeometry("geometry.root");

  // Load PWG1 library
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG1.so");
  
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  
  // Create the task
  AliAlignmentDataFilterITS *taskFilter = new AliAlignmentDataFilterITS("filterITS");
  AliLog::SetClassDebugLevel("AliAlignmentDataFilterITS",10);
  // configuration via AliITSRecoParam (should be taken from OCDB)
  AliITSReconstructor *itsReconstructor = new AliITSReconstructor();
  AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsReconstructor->SetRecoParam(itsRecoParam);


  // Add ESD handler
  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  TChain *chainESD = new TChain("esdTree");
  chainESD->Add("AliESDs.root");
 
  // Attach input
  cInput = mgr->CreateContainer("cInput",TChain::Class(),AliAnalysisManager::kInputContainer);
  //mgr->ConnectInput(taskFilter, 0, cInput); // v4-16-Release
  mgr->ConnectInput(taskFilter,0,mgr->GetCommonInputContainer());

  // Attach output
  cOutput0= mgr->CreateContainer("cOutput0",TTree::Class(),AliAnalysisManager::kOutputContainer,"AliTrackPoints.root");
  mgr->ConnectOutput(taskFilter,0,cOutput0);
  cOutput1= mgr->CreateContainer("cOutput1",TList::Class(),AliAnalysisManager::kOutputContainer,"AliTrackPoints.root");
  mgr->ConnectOutput(taskFilter,1,cOutput1);
  
  // Enable debug printouts
  mgr->SetDebugLevel(10);
  
  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();  
  mgr->StartAnalysis("local",chainESD,nentries,firstentry);

 return;
}

