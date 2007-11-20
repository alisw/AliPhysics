void runProtonAnalysis() {
  TStopwatch timer;
  timer.Start();
  
  //runLocal();
  //runInteractive();
  //runBatch();
  runProof();

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal() {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
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

  //____________________________________________________________//
  //_____________Setting up PWG2spectra.par_____________________//
  //____________________________________________________________//
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");

  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis("ESD"); 
  TagAna->ChainLocalTags("/home/pchrist/ALICE/Alien/Tutorial/November2007/Tags");

  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  
  TChain* chain = 0x0;
  chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *task1 = new AliAnalysisTaskProtons("TaskProtons");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Protons.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runInteractive() {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libProofPlayer.so");

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
 
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

  //____________________________________________________________//
  //_____________Setting up PWG2spectra.par_____________________//
  //____________________________________________________________//
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
  
  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
 
  //grid tags
  TAlienCollection* coll = TAlienCollection::Open("tag.xml");
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  TagAna->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *task1 = new AliAnalysisTaskProtons("TaskProtons");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Protons.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
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

  //____________________________________________________________//
  //_____________Setting up PWG2spectra.par_____________________//
  //____________________________________________________________//
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");

  //ANALYSIS PART
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
  const char *collectionfile = "wn.xml";

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *task1 = new AliAnalysisTaskProtons("TaskProtons");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Protons.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runProof() {
  TStopwatch timer;
  timer.Start();
  printf("****** Connect to PROOF *******\n");
  TProof::Open("proof://lxb6046.cern.ch"); 
  gProof->SetParallel(1);

  // Enable the Analysis Package
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("PWG2spectra.par");
  gProof->EnablePackage("PWG2spectra");
  
  // You should get this macro and the txt file from:
  // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
  gROOT->LoadMacro("CreateESDChain.C");
  TChain* chain = 0x0;
  chain = CreateESDChain("ESD82XX_30K.txt",10);
  chain->SetBranchStatus("*Calo*",0);

  gProof->Load("AliAnalysisTaskProtons.cxx++");

    //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *task1 = new AliAnalysisTaskProtons("TaskProtons");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Protons.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
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
