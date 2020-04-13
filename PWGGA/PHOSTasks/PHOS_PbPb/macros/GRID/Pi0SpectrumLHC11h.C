void Pi0SpectrumLHC11h(const char* dataset="collection.xml",
		       bool recompile = true
)
{
  /* $Id$ */
  TStopwatch timer;
  timer.Start();


    TStringToken libs("Core,Tree,Geom,VMC,Physics,Minuit,Gui,XMLParser,Minuit2,Proof,STEERBase,ESD,AOD,OADB,ANALYSIS,ANALYSISalice,CDB,RAWDatabase,STEER,CORRFW,PHOSUtils,PHOSbase,PHOSpi0Calib,PHOSrec,PHOSshuttle,PHOSsim", ",");
  while( libs.NextToken() )
    gSystem->Load( Form("lib%s", libs.Data()) );

  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  
  //load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");

  // A task can be compiled dynamically with AClic
  if( recompile ) {
    gROOT->LoadMacro("AliCaloPhoton.cxx+g");
    gROOT->LoadMacro("AliPHOSEPFlattener.cxx+g");
    gROOT->LoadMacro("AliAnalysisTaskPi0Flow.cxx+g");
  }
  else {
    gSystem->Load("libPWGGAPHOSTasks");
  }
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  if (1) // token == "OK" ) 
    TGrid::Connect("alien://");
  else 
    AliInfo("You are not connected to the GRID") ; 

  cout << "Pi0Analysis: processing collection " << dataset << endl;

  // Create the chain
  TChain* chain = new TChain("aodTree");

  TGridCollection * collection = gGrid->OpenCollection(dataset);
  
  TGridResult* result = collection->GetGridResult("", 0, 0);
  TList* rawFileList = result->GetFileInfoList();
  
  for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
    TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ; 
    const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;  
    printf("Processing %s\n", rawFile) ;
    chain->Add(rawFile);
    printf("Chain: %d entries.\n",chain->GetEntriesFast()); 
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  mgr->SetCommonFileName("histos.root");
  
  // AOD input handler
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  
  // Debug level
  mgr->SetDebugLevel(2);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
  AliEPSelectionTask *taskEP = AddTaskEventplane() ; 

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AliVZEROEPSelectionTask *selTask = AddTaskVZEROEPSelection();  

  // Add my task
  AliAnalysisTaskPi0Flow* task;
  if ( recompile ) {
    task = new AliAnalysisTaskPi0Flow("AliAnalysisTaskPi0Flow");
    //task->SetPeriod(AliAnalysisTaskPi0Flow::kLHC11h);
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("outCont1", TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root"));
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/AddTaskPHOSPi0Flow.C");
    task = AddTaskPHOSPi0Flow();
  }

/*
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
*/


  
  // // Create containers for input/output
  // AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  // AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("pi0SpecOut1",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  
  // // Connect input/output
  // mgr->ConnectInput(task1 , 0, cinput);
  // mgr->ConnectOutput(task1, 1, coutput1);
 
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }

  timer.Stop();
  timer.Print();
}
