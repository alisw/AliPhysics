void runLocal(const char *chainlistfile, int dataFromAlien=0) {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  if (dataFromAlien)
    TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");

  //____________________________________________________//
  //_____________Setting up STEERBase.par_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  //____________________________________________________//
  //_____________Setting up AOD.par_____________________//
  //____________________________________________________//
  setupPar("AOD");
  gSystem->Load("libAOD.so");

  //_________________________________________________________//
  //_____________Setting up ANALYSIS.par_____________________//
  //_________________________________________________________//
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");

  //_________________________________________________________//
  //_____________Setting up ANALYSISalice.par________________//
  //_________________________________________________________//
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

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
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I\"/usr/local/CERN/root/include\" -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I./ESD -I./AOD -I./ANALYSIS -I./PWG2AOD/AOD");
  gROOT->LoadMacro("ConfigFemtoAnalysis.C++");

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  //  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");

  if (dataFromAlien) {
    AliTagAnalysis *analysis = new AliTagAnalysis();
    chain = analysis->GetChainFromCollection(chainlistfile,"esdTree");
  }
  else {
    gROOT->LoadMacro("CreateESDChain.C");
    chain = CreateESDChain(chainlistfile,500);
  }

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  AliMCEventHandler *mcH = new AliMCEventHandler;

  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  
  mgr->SetMCtruthEventHandler(mcH);
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskFemto *task1 = new AliAnalysisTaskFemto("TaskFemto");

  mgr->AddTask(task1);

  // Create containers for input/output
  //  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("input0", 
							   TTree::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Femto.ESD.root");
  
  //____________________________________________//
  cinput1->SetData(chain);
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

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
