SteerAnalysisTaskTOFSpectraPbPb(const Char_t *inputfilename, Bool_t mcFlag = kFALSE, Bool_t mcTuneFlag = kFALSE, Bool_t pbpbFlag = kFALSE, Int_t maxFiles = kMaxInt, Int_t maxEv = kMaxInt)
{

  /* include path for ACLic */
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  /* load libraries */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  /* build analysis task class */
  gROOT->LoadMacro("AliAnalysisParticle.cxx+g");
  gROOT->LoadMacro("AliAnalysisEvent.cxx+g");
  gROOT->LoadMacro("AliAnalysisTrack.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskTOFSpectraPbPb.cxx+g");

  /* setup input chain */
  TString str = inputfilename;
  const Char_t *filename;
  TChain *chain = new TChain("esdTree");
  if (str.EndsWith(".xml")) {
    TGrid::Connect("alien://");
    Info("SteerTaskTOFSpectraPbPb", "reading data list from collection:");
    TGridCollection *coll = gGrid->OpenCollection(inputfilename, maxFiles);
    coll->Reset();
    while (coll->Next()) {
      filename = coll->GetTURL();
      Info("SteerTaskTOFSpectraPbPb", Form("%s", filename));
      chain->Add(filename);
    }
  }
  else if (str.EndsWith(".txt")) {
    Info("SteerTaskTOFSpectraPbPb", "reading data list from text file:");
    ifstream is(inputfilename);
    Char_t buf[4096];
    while(!is.eof()) {
      is.getline(buf, 4096);
      if (is.eof()) break;
      chain->Add(buf);
      Info("SteerTaskTOFSpectraPbPb", Form("%s", buf));
    }
    is.close();
  }
  else {
    Info("SteerTaskTOFSpectraPbPb", "single file:");
    filename = inputfilename;
    Info("SteerTaskTOFSpectraPbPb", Form("%s", filename));
    chain->Add(filename);
  }
  Info("SteerTaskTOFSpectraPbPb", Form("chain is ready: %d events", chain->GetEntries()));

  /* create analysis manager */
  AliAnalysisManager *mgr = new AliAnalysisManager("TOFSpectraPbPb");

  /* define input event handler */
  AliESDInputHandler *esdh = new AliESDInputHandler();
  esdh->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdh);

  /* define MC truth event handler */
  if (mcFlag) {
    AliMCEventHandler *mch = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mch);
  }

  /* define output handler */
  AliAODHandler *outputh = new AliAODHandler();
  mgr->SetOutputEventHandler(outputh);

  /* add tasks */
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(mcFlag);
  if (pbpbFlag) {
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality(); 
    //    centralityTask->SetPass(2);
    if (mcFlag) centralityTask->SetMCInput();
  }
  gROOT->LoadMacro("AddAnalysisTaskTOFSpectraPbPb.C");
  AliAnalysisTaskTOFSpectraPbPb *thisTask = AddAnalysisTaskTOFSpectraPbPb(mcFlag, mcTuneFlag, pbpbFlag);

  /* start analysis */
  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, maxEv);

  /* create dummy file to tell we are done */
  gSystem->Exec("touch done");

}
