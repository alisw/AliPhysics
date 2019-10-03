SteerAnalysisTaskParticleEfficiencyLF(const Char_t *inputfilename, Int_t maxFiles = kMaxInt, Int_t maxEv = kMaxInt)
{

  /* include path for ACLic */
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  /* load libraries */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  /* build analysis task class */
  gROOT->LoadMacro("AliAnalysisTaskParticleEfficiencyLF.cxx+g");

  /* setup input chain */
  TString str = inputfilename;
  const Char_t *filename;
  TChain *chain = new TChain("esdTree");
  if (str.EndsWith(".xml")) {
    TGrid::Connect("alien://");
    Info("SteerTaskParticleEfficiency", "reading data list from collection:");
    TGridCollection *coll = gGrid->OpenCollection(inputfilename, maxFiles);
    coll->Reset();
    while (coll->Next()) {
      filename = coll->GetTURL();
      Info("SteerTaskParticleEfficiency", Form("%s", filename));
      chain->Add(filename);
    }
  }
  else if (str.EndsWith(".txt")) {
    Info("SteerTaskParticleEfficiency", "reading data list from text file:");
    ifstream is(inputfilename);
    Char_t buf[4096];
    while(!is.eof()) {
      is.getline(buf, 4096);
      if (is.eof()) break;
      chain->Add(buf);
      Info("SteerTaskParticleEfficiency", Form("%s", buf));
    }
    is.close();
  }
  else {
    Info("SteerTaskParticleEfficiency", "single file:");
    filename = inputfilename;
    Info("SteerTaskParticleEfficiency", Form("%s", filename));
    chain->Add(filename);
  }
  Info("SteerTaskParticleEfficiency", Form("chain is ready: %d events", chain->GetEntries()));

  /* create analysis manager */
  AliAnalysisManager *mgr = new AliAnalysisManager("ParticleEfficiency");

  /* define input event handler */
  AliESDInputHandler *esdh = new AliESDInputHandler();
  esdh->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdh);

  /* define MC truth event handler */
  AliMCEventHandler *mch = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mch);
  
  /* define output handler */
  AliAODHandler *outputh = new AliAODHandler();
  mgr->SetOutputEventHandler(outputh);
  
  /* add tasks */
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE);
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *centralityTask = AddTaskCentrality(); 
  centralityTask->SetMCInput();
  gROOT->LoadMacro("AddAnalysisTaskParticleEfficiency.C");
  
  AddAnalysisTaskParticleEfficiency("pi+");
  AddAnalysisTaskParticleEfficiency("pi-");
  AddAnalysisTaskParticleEfficiency("K+");
  AddAnalysisTaskParticleEfficiency("K-");
  AddAnalysisTaskParticleEfficiency("proton");
  AddAnalysisTaskParticleEfficiency("antiproton");
  AddAnalysisTaskParticleEfficiency("K*0");
  AddAnalysisTaskParticleEfficiency("K*0_bar");
  AddAnalysisTaskParticleEfficiency("phi");

  /* start analysis */
  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, maxEv);

  /* create dummy file to tell we are done */
  gSystem->Exec("touch done");

}
