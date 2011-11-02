void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  // Use par files only for PWG2 code
  int useParFiles = 0;
  int usePWG2ParFiles = 1;

  // Use precompiled libraries for the analysis framework
  if (useParFiles)
    setupPar("STEERBase");
  if (gSystem->Load("libSTEERBase.so")<0) {
    cout << "Cannot load libSTEERBase.so . Exiting" << endl;
    exit(0);
  }
  gSystem->Load("libVMC.so");

  if (useParFiles)
    setupPar("ESD");
  if (gSystem->Load("libESD.so")<0) {
    cout << "Cannot load libESD.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("AOD");
  if (gSystem->Load("libAOD.so")<0) {
    cout << "Cannot load libAOD.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("ANALYSIS");
  if (gSystem->Load("libANALYSIS.so")<0) {
    cout << "Cannot load libANALYSIS.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("ANALYSISalice");
  if (gSystem->Load("libANALYSISalice.so")<0) {
    cout << "Cannot load libANALYSISalice.so . Exiting" << endl;
    exit(0);
  }

  //____________________________________________________//
  //_____________Setting up PWG2AOD.par_________________//
  //____________________________________________________//
  if (usePWG2ParFiles) {
    // char dynpath[10000];
    // sprintf(dynpath, ".:%s", gSystem->GetDynamicPath());
    // gSystem->SetDynamicPath(dynpath);
    TString dynpath;
    dynpath = ".:";
    dynpath += gSystem->GetDynamicPath();
    gSystem->SetDynamicPath(dynpath.Data());
  }

  if (usePWG2ParFiles) {
    setupPar("PWG2AOD");
    if (gSystem->Load("./PWG2AOD/libPWG2AOD.so")<0) {
      cout << "Cannot load local libPWG2AOD.so . Exiting" << endl;
      exit(0);
    }
  }
  else {
    if (gSystem->Load("libPWG2AOD.so")<0) {
      cout << "Cannot load libPWG2AOD.so . Exiting" << endl;
      exit(0);
    }
  }
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopy.par__________//
  //____________________________________________________//
  if (usePWG2ParFiles) {
    setupPar("PWG2femtoscopy");
    if (gSystem->Load("./PWG2femtoscopy/libPWG2femtoscopy.so")<0) {
      cout << "Cannot load local libPWG2femtoscopy.so . Exiting" << endl;
      exit(0);
    }
  }
  else {
    if (gSystem->Load("libPWG2femtoscopy.so")<0) {
      cout << "Cannot load libPWG2femtoscopy.so . Exiting" << endl;
      exit(0);
    }
  }
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopyUser.par______//
  //____________________________________________________//
  if (usePWG2ParFiles) {
    setupPar("PWG2femtoscopyUser");
    if (gSystem->Load("./PWG2femtoscopyUser/libPWG2femtoscopyUser.so")<0) {
      cout << "Cannot load libPWG2femtoscopyUser.so . Exiting" << endl;
      exit(0);
    }
  }
  else {
    if (gSystem->Load("libPWG2femtoscopyUser.so")<0) {
      cout << "Cannot load libPWG2femtoscopyUser.so . Exiting" << endl;
      exit(0);
    }
  }
  
  //ANALYSIS PART
  const char *collectionfile="wn.xml";
  
  //____________________________________________//
  //Usage of event tags
  // AliTagAnalysis *analysis = new AliTagAnalysis();
  // TChain *chain = 0x0;
  // chain = analysis->CreateChainFromCollection(collectionfile,"esdTree");

  TChain *chain = new TChain("aodTree");;
  //  gROOT->LoadMacro("CreateESDChain.C");
//  const char *collectionfile="/home/akisiel/LHC10h.esds.txt";
  
  ifstream *istr = new ifstream(collectionfile);
  
  char fname[2000];
  char pname[2000];
  while (!istr->eof()) {
    fname[0] = '\0';
    (*istr) >> fname;
    if (strlen(fname) > 10) {
      sprintf(pname, "alien://%s", fname);
      chain->Add(pname);
      
    }
  }
  
//  chain->Add("data/AliAOD.root");

  //  chain->Add("Data121040/AliESDs.root");
  
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/001/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/002/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/003/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/004/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/005/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/006/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/007/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/008/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/009/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/010/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/011/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/012/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/013/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/014/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/015/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/016/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/017/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/018/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/019/AliESDs.root");
//   chain->Add("alien:///alice/sim/LHC10g1a/130844/020/AliESDs.root");
//  chain->Add("ESDs/AliESDs.101.root");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
//   AliESDInputHandler* esdH = new AliESDInputHandler;
//  esdH->SetInactiveBranches("FMD CaloCluster");
//   esdH->SetReadFriends(kFALSE);
  AliAODInputHandler *aodH = new AliAODInputHandler;
  //  mgr->SetInputEventHandler(esdH);  
  mgr->SetInputEventHandler(aodH);

//  AliMCEventHandler *mcH = new AliMCEventHandler;
//  mgr->SetMCtruthEventHandler(mcH);

//  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
//  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(0,  0);
  //physSelTask->GetPhysicsSelection()->SetUseBXNumbers(kFALSE);

  // AOD output handler
//   AliAODHandler* aodHandler   = new AliAODHandler();
//   aodHandler->SetOutputFileName("AliAOD.root");
//   //aodHandler->SetFillAOD(kFALSE);
//   aodHandler->SetFillAODforRun(kFALSE);
//   mgr->SetOutputEventHandler(aodHandler);
  
//   AliCentralitySelectionTask *centralityTask = new AliCentralitySelectionTask("CentralitySelection");
//   centralityTask->SetPass(2);
//   mgr->AddTask(centralityTask);
//   mgr->ConnectInput(centralityTask, 0, mgr->GetCommonInputContainer());
  
//   AliAnalysisDataContainer *outputCentrality =  mgr->CreateContainer("outputCentrality", TList::Class(),
//   AliAnalysisManager::kOutputContainer, "CentralityOutput.root");
//   mgr->ConnectOutput(centralityTask, 1, outputCentrality);

  //  ESD filter task configuration.
  //  gROOT->LoadMacro("AddTaskESDFilter.C");
  // TODO usually requires muon filter libs
  //  AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kFALSE,  kFALSE, kFALSE, kTRUE, kTRUE);

  //AddTaskPIDResponse
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *pidresponse = AddTaskPIDResponse();

  //____________________________________________//
  // 1st Pt task
  gROOT->LoadMacro("AddTaskFemto.C");
  AliAnalysisTaskFemto *taskfemto = AddTaskFemto("./ConfigFemtoAnalysis.C");

  //____________________________________________//
  // Run the analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runProcess","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }

  return 1;
}
