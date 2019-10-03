SteerAnalysisTaskPIDFluctuation(const Char_t *inputfilename, Int_t maxFiles = kMaxInt, Int_t maxEv = kMaxInt)
{

  /* include path for ACLic */
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  /* load libraries */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  /* build analysis task class */
  gROOT->LoadMacro("AliAnalysisTaskPIDFluctuation.cxx+g");

  /* setup input chain */
  TString str = inputfilename;
  const Char_t *filename;
  TChain *chain = new TChain("esdTree");
  if (str.EndsWith(".xml")) {
    TGrid::Connect("alien://");
    Info("", "reading data list from collection:");
    TGridCollection *coll = gGrid->OpenCollection(inputfilename, maxFiles);
    coll->Reset();
    while (coll->Next()) {
      filename = coll->GetTURL();
      Info("", Form("%s", filename));
      chain->Add(filename);
    }
  }
  else if (str.EndsWith(".txt")) {
    Info("", "reading data list from text file:");
    ifstream is(inputfilename);
    Char_t buf[4096];
    while(!is.eof()) {
      is.getline(buf, 4096);
      if (is.eof()) break;
      chain->Add(buf);
      Info("", Form("%s", buf));
    }
    is.close();
  }
  else {
    Info("", "single file:");
    filename = inputfilename;
    Info("", Form("%s", filename));
    chain->Add(filename);
  }
  Info("", Form("chain is ready: %d events", chain->GetEntries()));

  /* create analysis manager */
  AliAnalysisManager *mgr = new AliAnalysisManager("EbyEFluctuationPID");

  /* define input event handler */
  AliESDInputHandler *esdh = new AliESDInputHandler();
  esdh->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdh);

  /* add tasks */
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE);
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *centralityTask = AddTaskCentrality(); 
  gROOT->LoadMacro("AddAnalysisTaskPIDFluctuation.C");
  AliAnalysisTaskPIDFluctuation *thisTask = AddAnalysisTaskPIDFluctuation();

  printf("ready-steady-go\n");

  /* start analysis */
  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, maxEv);

  /* create dummy file to tell we are done */
  gSystem->Exec("touch done");

}
