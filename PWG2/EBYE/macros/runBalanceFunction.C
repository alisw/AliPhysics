void runBalanceFunction() {
  TStopwatch timer;
  timer.Start();

  runLocal();

  timer.Stop();
  timer.Print();
}

void runLocal() {
  //____________________________________________________//
  //_____________Setting up the par files_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  setupPar("AOD");
  gSystem->Load("libAOD.so");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");
  setupPar("PWG2ebye");
  gSystem->Load("libPWG2ebye.so");

  //___________________________________________________//
  TChain *chain = new TChain("esdTree");
  chain->Add("AliESDs.root");


  //___________________________________________________//
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Create task
  gROOT->LoadMacro("AliAnalysisTaskBF.cxx+g");
  AliAnalysisTaskBF *task = new AliAnalysisTaskBF("TaskBF");

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("contBF", AliBalance::Class(), AliAnalysisManager::kOutputContainer, "BF.ESD.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}

//___________________________________________//
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
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
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
