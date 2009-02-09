//--------------------------------------------------------------------------
// Base macro for submitting AliMUONRecoCheck analysis.
//
// The macro reads ESDs, Kinematics and TrackRefs and outputs file:
// - RecoCheck.root
//--------------------------------------------------------------------------

void RunRecoCheck(Bool_t local = kFALSE) {

  TStopwatch timer;
  timer.Start();

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  gSystem->Load("libPWG3muondep.so");

  TChain* chain = new TChain("esdTree");

  if (!local) {

    printf("*** Connect to AliEn ***\n");
    TGrid::Connect("alien://");
    gSystem->Load("libProofPlayer.so");

    TAlienCollection* coll = TAlienCollection::Open("wn.xml");
    TGridResult* result = coll->GetGridResult("",0,0);
    for(Int_t i = 0; i < result->GetEntries(); i++) {
      printf("TURL = %s \n",result->GetKey(i,"turl"));
      chain->Add(result->GetKey(i,"turl"));
    }
  } else {
    chain->Add("/path_1_to/AliESDs.root");
    chain->Add("/path_2_to/AliESDs.root");
    chain->Add("/path_3_to/AliESDs.root");
  }

  // set the magnetic field for track extrapolations
  // look in the OCDB for the value of the magnet current
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //cdb->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-10-Release/Ideal/");
  cdb->SetRun(0);

  AliCDBEntry *entry;
  entry = cdb->Get("GRP/GRP/Data");
  AliGRPObject *obj = (AliGRPObject*)entry->GetObject();
  Float_t currentL3 = obj->GetL3Current(0);
  printf("Loading field map for L3 current = %f A ...\n",currentL3);

  AliLog::SetGlobalLogLevel(AliLog::kError);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);

  //____________________________________________//
  // ntuple task
  AliAnalysisTaskRecoCheck *task = new AliAnalysisTaskRecoCheck("TaskRecoCheck");
  task->SetL3Current(currentL3);
  mgr->AddTask(task);
  
  // Create containers for input/output

  // input
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain",TChain::Class(),AliAnalysisManager::kInputContainer);

  // output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("ctree", TTree::Class(),AliAnalysisManager::kOutputContainer,"RecoCheck.root");

  //____________________________________________//
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,0,coutput);

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();

    mgr->StartAnalysis("local",chain);

  }

  timer.Stop();
  timer.Print();

}

