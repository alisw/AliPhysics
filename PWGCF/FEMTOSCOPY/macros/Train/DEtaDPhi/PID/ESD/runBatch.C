void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  //TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  // Use par files only for PWG2 code
  int useParFiles = 0;
  int usePWGCFParFiles = 0;

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
  if (usePWGCFParFiles) {
    // char dynpath[10000];
    // sprintf(dynpath, ".:%s", gSystem->GetDynamicPath());
    // gSystem->SetDynamicPath(dynpath);
        TString dynpath;
    dynpath = ".:";
    dynpath += gSystem->GetDynamicPath();
    gSystem->SetDynamicPath(dynpath.Data());
    }

    //if (usePWGCFParFiles) {
    // setupPar("PWGCFAOD");
    //if (gSystem->Load("./PWGCFAOD/libPWGCFAOD.so")<0) {
    //  cout << "Cannot load local libPWGCFAOD.so . Exiting" << endl;
    //  exit(0);
    // }
    // }
    //else {
    //if (gSystem->Load("libPWGCFAOD.so")<0) {
    //  cout << "Cannot load libPWGCFAOD.so . Exiting" << endl;
        //exit(0);
        //}
        //}
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopy.par__________//
  //____________________________________________________//
  if (usePWGCFParFiles) {
    setupPar("PWGCFfemtoscopy");
    if (gSystem->Load("./PWGCFfemtoscopy/libPWGCFfemtoscopy.so")<0) {
      cout << "Cannot load local libPWG2femtoscopy.so . Exiting" << endl;
      exit(0);
    }
  }
  else {
    if (gSystem->Load("libPWGCFfemtoscopy.so")<0) {
      cout << "Cannot load libPWGCFfemtoscopy.so . Exiting" << endl;
      exit(0);
    }
  }
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopyUser.par______//
  //____________________________________________________//
  if (usePWGCFParFiles) {
    setupPar("PWGCFfemtoscopyUser");
    if (gSystem->Load("./PWGCFfemtoscopyUser/libPWGCFfemtoscopyUser.so")<0) {
      cout << "Cannot load libPWGCFfemtoscopyUser.so . Exiting" << endl;
      exit(0);
    }
  }
  else {
    if (gSystem->Load("libPWGCFfemtoscopyUser.so")<0) {
      cout << "Cannot load libPWGCFfemtoscopyUser.so . Exiting" << endl;
      exit(0);
    }
  }
  
  //ANALYSIS PART
  const char *collectionfile="wn.xml";
  
  //____________________________________________//
  //Usage of event tags


    TChain *chain = new TChain("esdTree");
    //chain->Add("../Data/PbPb/AOD/1/AliAOD.root");

    //gROOT->LoadMacro("./CreateESDChain.C");
    //TChain *chain = CreateESDChain("files.txt");
  
  // ifstream *istr = new ifstream(collectionfile);
  
  // char fname[2000];
  // char pname[2000];
  // while (!istr->eof()) {
  //   fname[0] = '\0';
  //   (*istr) >> fname;
  //   if (strlen(fname) > 10) {
  //     sprintf(pname, "alien://%s", fname);
  //     chain->Add(pname);
      
  //   }
  //   }

  chain->Add("/opt/alice/workdir/TestConfig/data/ESD/pp/LHC10e/pass2/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/1/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/2/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/3/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/4/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/5/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/6/AliESDs.root");
    // chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/7/AliESDs.root");
    //chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/8/AliESDs.root");
    //chain->Add("/opt/alice/workdir/TestConfig/data/Pythia/9/AliESDs.root");
  //____________________________________________//
  // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
    //AliESDInputHandler* esdH = new AliESDInputHandler;
//  esdH->SetInactiveBranches("FMD CaloCluster");
    
    AliESDInputHandler *esdH = new AliESDInputHandler;
    esdH->SetReadFriends(kFALSE);
  //  mgr->SetInputEventHandler(esdH);  
  mgr->SetInputEventHandler(esdH);

//  AliMCEventHandler *mcH = new AliMCEventHandler;
//  mgr->SetMCtruthEventHandler(mcH);


  //  ESD filter task configuration.
  //  gROOT->LoadMacro("AddTaskESDFilter.C");
  // TODO usually requires muon filter libs
  //  AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kFALSE,  kFALSE, kFALSE, kTRUE, kTRUE);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();

  //AddTaskPIDResponse
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *pidresponse = AddTaskPIDResponse();
 
  //____________________________________________//
  // 1st Pt task
  gROOT->LoadMacro("AddTaskFemto.C");
  AliAnalysisTaskFemto *taskfemto = AddTaskFemto("./ConfigFemtoAnalysis.C");
  //taskfemto->SelectCollisionCandidates(AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kMB);
  taskfemto->SelectCollisionCandidates(AliVEvent::kMB);
    
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
