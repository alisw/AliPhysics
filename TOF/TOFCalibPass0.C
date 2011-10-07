//___________________________________________________________________________

LoadLibraries()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTOFcalib");
}

//___________________________________________________________________________

MakeOCDB(const Char_t *filename = "TOFCalibPass0.root", const Char_t *dbString = "local://$HOME/OCDB")
{
  LoadLibraries();
  AliTOFAnalysisTaskCalibPass0 calibTask;
  calibTask.ProcessOutput(filename, dbString);
  printf("TOF calibration status code: %d\n", calibTask.GetStatus()); 
}

//___________________________________________________________________________

SteerTask(const Char_t *inputfilename, Int_t maxFiles = kMaxInt, Int_t maxEv = kMaxInt)
{

  LoadLibraries();

  /* setup input chain */
  TString str = inputfilename;
  const Char_t *filename;
  TChain *chain = new TChain("esdTree");
  if (str.EndsWith(".xml")) {
    TGrid::Connect("alien://");
    Info("SteerTaskEventTime", "reading data list from collection:");
    TAlienCollection coll(inputfilename, maxFiles);
    coll.Reset();
    while (coll.Next()) {
      filename = coll.GetTURL();
      Info("SteerTaskEventTime", Form("%s", filename));
      chain->Add(filename);
    }
  }
  else if (str.EndsWith(".txt")) {
    Info("SteerTaskEventTime", "reading data list from text file:");
    ifstream is(inputfilename);
    Char_t buf[4096];
    while(!is.eof()) {
      is.getline(buf, 4096);
      if (is.eof()) break;
      chain->Add(buf);
      Info("SteerTaskEventTime", Form("%s", buf));
    }
    is.close();
  }
  else {
    Info("SteerTaskEventTime", "single file:");
    filename = inputfilename;
    Info("SteerTaskEventTime", Form("%s", filename));
    chain->Add(filename);
  }
  Info("SteerTaskEventTime", Form("chain is ready: %d events", chain->GetEntries()));

  /* create analysis manager */
  AliAnalysisManager *mgr = new AliAnalysisManager("EventTime");

  /* define input event handler */
  AliESDInputHandler *esdh = new AliESDInputHandler();
  esdh->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdh);

  /* add tasks */
  gROOT->LoadMacro("$ALICE_ROOT/TOF/AddTOFAnalysisTaskCalibPass0.C");
  AliTOFAnalysisTaskCalibPass0 *thisTask = AddTOFAnalysisTaskCalibPass0();

  /* start analysis */
  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, maxEv);

  /* create dummy file to tell we are done */
  gSystem->Exec("touch done");

}
