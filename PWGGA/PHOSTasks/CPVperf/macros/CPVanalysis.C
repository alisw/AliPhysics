void CPVanalysis(const TString dataset="")
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
  gROOT->LoadMacro("AliAnalysisTaskCPV.cxx+g");
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  token="OK";
  if ( token == "OK" ) 
    TGrid::Connect("alien://");
  else 
    Printf("You are not connected to the GRID") ; 
    // AliInfo("You are not connected to the GRID") ; 

  cout << "CPVanalysis: processing collection " << dataset << endl;
  
  // Create the chain
    TChain* chain;

  if (dataset.Contains(".txt")) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    chain = CreateESDChain(dataset.Data(), 1000);
  }
  else if (dataset.Contains(".xml")) {
    chain = new TChain("esdTree");
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
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("CPVanalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // Debug level
  mgr->SetDebugLevel(0);
  
  //PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTask *taskPID =  AddTaskPIDResponse(/*Bool_t isMC=*/ kFALSE, /*Bool_t autoMCesd=*/kTRUE,
						 /*Bool_t tuneOnData=*/kTRUE);

  // // Add physics selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  physSelTask->GetPhysicsSelection()->SetUseBXNumbers(kFALSE);// ---> ???

  // Add my task
  AliAnalysisTaskCPV *task1 = new AliAnalysisTaskCPV("CPVanalysis");

  // task1->SelectCollisionCandidates(AliVEvent::kMB);
  mgr->AddTask(task1);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histCPV1",TList::Class(),AliAnalysisManager::kOutputContainer,"CPVanalysis.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histCPV2",TList::Class(),AliAnalysisManager::kOutputContainer,"CPVanalysis.root");

  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(1);


  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  mgr->ConnectOutput(task1, 2, coutput2);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
}
