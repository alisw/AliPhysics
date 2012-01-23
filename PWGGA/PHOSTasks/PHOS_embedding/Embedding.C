void Embedding(const char* dataset="collection.xml")
{
    
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");

  // A task can be compiled dynamically with AClic
  gROOT->LoadMacro("AliPHOSEmbedding.cxx+g");
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  if (1) // token == "OK" ) 
    TGrid::Connect("alien://");
  else 
    AliInfo("You are not connected to the GRID") ; 
  cout << "Pi0Analysis: processing collection " << dataset << endl;
  
  // Create the chain
  TChain* chain = new TChain("esdTree");

  chain->Add("/home/prsnko/PbPb/Embedding/Data/AliESDs.root") ;
/*

  TGridCollection * collection = dynamic_cast<TGridCollection*>(TAlienCollection::Open(dataset));
  
  TAlienResult* result = collection->GetGridResult("",0 ,0);
  TList* rawFileList = result->GetFileInfoList();
  for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
    TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ; 
    const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;  
    printf("Processing %s\n", rawFile) ;
    chain->Add(rawFile);
    printf("Chain: %d entries.\n",chain->GetEntries()); 
  }
  TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(0));
  const char * fn = fi->GetCurrentUrl()->GetUrl() ;
*/
  char runNum[7]={"139438"}; 
//  for(Int_t i=0;i<6;i++)runNum[i]=fn[35+i] ;
  runNum[6]=0 ;
  Int_t iRunNum=atoi(runNum) ;
  printf("Run number=%d \n",iRunNum) ;

  //Run AOD simulation
  int nrun = atoi(runNum);
  int nevent = 0;
  int seed = 0;

  char sseed[1024];
  char sevent[1024];
  char sprocess[1024];
  char sfield[1024];
  char senergy[1024];

  sprintf(sevent,"");
  sprintf(sprocess,"");
  sprintf(sfield,"");
  sprintf(senergy,"");

  seed = 0;
  sprintf(sseed,"%d",seed);

  if (seed==0) {
    fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr,"!!!!  WARNING! Seeding variable for MC is 0          !!!!\n");
    fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  } else {
    fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stdout,"!!!  MC Seed is %d \n",seed);
    fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
  
// set the seed environment variable
  gSystem->Setenv("CONFIG_SEED",sseed);
  gSystem->Setenv("CONFIG_RUN_TYPE","kPythia6"); // kPythia6 or kPhojet^M
  gSystem->Setenv("CONFIG_FIELD","k5kG");      // kNoField or k5kG^M
  gSystem->Setenv("CONFIG_ENERGY","2760");    // 900 or 10000 (GeV)
  gSystem->Setenv("DC_RUN",runNum);    //run number 
  

  char nSimEvents[55] ;
  sprintf(nSimEvents,"%d",chain->GetEntries());
  gSystem->Setenv("SIM_EVENTS",nSimEvents); 
  gSystem->Exec("mv geometry.root geometry_PHOS.root") ;
//  gSystem->Exec("aliroot -b -q simrun.C > simrun.log 2>&1");
  gSystem->Exec("mv geometry_PHOS.root geometry.root") ;

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0EmbeddingManager");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  // Output
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAODout.root");
  mgr->SetOutputEventHandler(aodHandler);

  
  // Debug level
  mgr->SetDebugLevel(0);

  // Add physics selection
//  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
//  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE,kFALSE);

  //Add centrality
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality() ;
  taskCentrality->SetPass(2); // remember to set the pass you are processing!!! 

  // Add my task
  AliPHOSEmbedding *task1 = new AliPHOSEmbedding("Embedding");

  TChain* chainAOD = new TChain("aodTree");
  chainAOD->AddFile("AliAOD.root") ;
  task1->SetSignalChain(chainAOD) ;
 // task1->SelectCollisionCandidates();


  TFile *fOldCalib = TFile::Open("OldCalibration.root");
  if(fOldCalib->IsOpen()){
    printf("\n\n...Adding PHOS calibration used in ESD production \n") ;
    char key[55] ;
    TH2F * hCalib[5] ;
    for(Int_t mod=0;mod<5; mod++){
      sprintf(key,"calibrationMod%d",mod) ;
      hCalib[mod] = (TH2F*)fOldCalib->Get(key) ;
    }
    task1->SetOldCalibration(hCalib) ;
    fOldCalib->Close() ;
  }

  mgr->AddTask(task1);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("output0",
                      TTree::Class(), AliAnalysisManager::kOutputContainer);


  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 0,coutput1);

  AliLog::SetClassDebugLevel("AliGeomManager", 10) ;

  AliCDBManager::Instance()->SetRun(iRunNum) ;
//  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB") ;
//  AliCDBManager::Instance()->SetDefaultStorage("raw://") ;
  AliCDBManager::Instance()->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB") ;
  AliCDBManager::Instance()->SetSpecificStorage("PHOS/*/*",
                                                "local://OCDB");
printf("RunNunm===%d \n",iRunNum) ;
  

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
  if(iRunNum<=137848){ //period 1
    gSystem->Exec("mv BadMap_LHC10h_period1.root BadMap_LHC10h.root") ;
  }
  else{
    gSystem->Exec("mv BadMap_LHC10h_period234.root BadMap_LHC10h.root") ;
  }
//  gSystem->Exec("aliroot -b -q Analyze.C> analyze.log 2>&1");

}
