void Pi0Spectrum(const char* dataset="")
{
    
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  
  //load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");

  // A task can be compiled dynamically with AClic
  gROOT->LoadMacro("AliCaloPhoton.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskPi0.cxx+g");
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  if ( token == "OK" ) 
    TGrid::Connect("alien://");
  else 
    AliInfo("You are not connected to the GRID") ; 

  cout << "Pi0Analysis: processing collection " << dataset << endl;
  
  // Create the chain
  TChain* chain = new TChain("esdTree");
  TGridCollection * collection = gGrid->OpenCollection(dataset);
  
  TGridResult* result = collection->GetGridResult("", 0, 0);
  TList* rawFileList = result->GetFileInfoList();
  
  for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
    TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ; 
    const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;  
    printf("Processing %s\n", rawFile) ;
    chain->Add(rawFile);
    printf("Chain: %d entries.\n",chain->GetEntries()); 
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // Debug level
  mgr->SetDebugLevel(0);

  // Add physics selection
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

  // Add my task
  AliAnalysisTaskPi0 *task1 = new AliAnalysisTaskPi0("Pi0Spectrum");

  TFile *fBadMap = TFile::Open("BadMap_LHC11a_pp2760_20130902.root");
  if(fBadMap->IsOpen()){
    printf("\n\n...Adding PHOS bad channel map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fBadMap->Get(key) ;
      if(h)
	task1->SetPHOSBadMap(mod,h) ;
    }
    fBadMap->Close() ;
  }
  // Set BC gap for LHC10e in seconds
  task1->SetBCgap(525.e-09);

  // Set abs.recalibration for LHC11a
  task1->SetRecalib(1, 0.9942);
  task1->SetRecalib(2, 0.9822);
  task1->SetRecalib(3, 1.0072);

  task1->SelectCollisionCandidates(AliVEvent::kMB);
  mgr->AddTask(task1);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histESD",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  
  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
}
