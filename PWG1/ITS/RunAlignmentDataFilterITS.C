void RunAlignmentDataFilterITS() {
  //
  // Macro to extract AliTrackPoints for ITS
  // A.Dainese, andrea.dainese@pd.infn.it
  //

  // Input
  Bool_t singlefile=kFALSE;
  //TString esdpath="/home/dainesea/alignData/RAWdata_CosmicsSum09/RecoSPD/chunk.";
  TString esdpath="/home/dainesea/alignData/RAWdata_CosmicsSum09/RecoITS_B_mille_SPD_SDDSSDsurvey_SSDHLayer_th50_130709/chunk.";
  Int_t ifirst=1, ilast=11;
  //
  Int_t nentries=1234567890;
  Int_t firstentry=0;


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
  AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetAlignFilterUseLayer(0,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(1,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(2,kFALSE);
  itsRecoParam->SetAlignFilterUseLayer(3,kFALSE);
  itsRecoParam->SetAlignFilterUseLayer(4,kFALSE);
  itsRecoParam->SetAlignFilterUseLayer(5,kFALSE);
  taskFilter->SetITSRecoParam(itsRecoParam);
  taskFilter->SetOnlySPDFO();

  // Add ESD handler
  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  TChain *chainESD = new TChain("esdTree");
  if(singlefile) {
    chainESD->Add("AliESDs.root");
  } else {
    for(Int_t i=ifirst; i<=ilast; i++) {
      TString esdfile=esdpath; esdfile+=i; esdfile.Append("/AliESDs.root");
      chainESD->Add(esdfile.Data());
    }
  } 

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
