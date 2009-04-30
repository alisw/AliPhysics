void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");

  // Use precompiled libraries for the analysis framework
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");

  // Use par files only for PWG2 code

  //____________________________________________________//
  //_____________Setting up PWG2AOD.par_________________//
  //____________________________________________________//
  setupPar("PWG2AOD");
  gSystem->Load("libPWG2AOD.so");
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopy.par__________//
  //____________________________________________________//
  setupPar("PWG2femtoscopy");
  gSystem->Load("libPWG2femtoscopy.so");
  
  //____________________________________________________//
  //_____________Setting up PWG2femtoscopyUser.par______//
  //____________________________________________________//
  setupPar("PWG2femtoscopyUser");
  gSystem->Load("libPWG2femtoscopyUser.so");
  
  //ANALYSIS PART
  const char *collectionfile="wn.xml";

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  

  AliMCEventHandler *mcH = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(mcH);

  //____________________________________________//
  // 1st Pt task
  gROOT->LoadMacro("AddTaskFemto.C");
  AliAnalysisTaskFemto *taskfemto = AddTaskFemto();

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
