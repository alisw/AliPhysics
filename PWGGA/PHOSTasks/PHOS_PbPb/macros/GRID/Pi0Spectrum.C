void Pi0Spectrum(const char* dataset="collection.xml")
{

  /* $Id$ */
    
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
  gROOT->LoadMacro("AliPHOSEPFlattener.cxx++") ;
  gROOT->LoadMacro("AliAnalysisTaskPi0Flow.cxx+g");
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  if (1) // token == "OK" ) 
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

  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality() ;
 // taskCentrality->SetPass(2); // remember to set the pass you are processing!!! 

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
  AliEPSelectionTask *taskEP = AddTaskEventplane() ; 
  
  // Add my task
  AliAnalysisTaskPi0Flow *task1 = new AliAnalysisTaskPi0Flow("Pi0Spectrum");

  TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10h_period1.root");
  //  TFile *fBadMap = TFile::Open("BadMap_LHC10h_period1.root");
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

  task1->SelectCollisionCandidates();
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
