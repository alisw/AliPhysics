void runProtonAnalysisQA() {
  TStopwatch timer;
  timer.Start();
  
  runProof(200000,"/COMMON/COMMON/LHC08c11_10TeV_0.5T"); //use data sets
  //runInteractive("wn.xml");
  
  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runInteractive(const char *collectionfile) {
  TString outputFilename1 = "Protons.QA.root"; 
  TString outputFilename2 = "Protons.MC.QA.root"; 
  TString outputFilename3 = "Protons.QA.Histograms.root"; 
  TString outputFilename4 = "Protons.Efficiency.root"; 

  TGrid::Connect("alien://");

  //Setup the par files
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");
  setupPar("ESD");
  gSystem->Load("libESD.so");
  setupPar("AOD");
  gSystem->Load("libAOD.so");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");
  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");

  gROOT->LoadMacro("AliAnalysisTaskProtonsQA.cxx+");
  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtonsQA *taskProtonsQA = new AliAnalysisTaskProtonsQA("TaskProtonsQA");
  taskProtonsQA->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
  taskProtonsQA->SetAcceptedVertexDiamond(5.,5.,15.);
  mgr->AddTask(taskProtonsQA);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
							   TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("globalQAList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename1.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("pdgCodeList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename2.Data());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("mcProcessList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename2.Data());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("acceptedCutList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename3.Data());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("acceptedDCAList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename3.Data());
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("efficiencyList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename4.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}
 
//_________________________________________________//
void runProof(Int_t stats = 0, const char* dataset = 0x0) {
  TString outputFilename1 = "Protons.QA.root"; 
  TString outputFilename2 = "Protons.MC.QA.root"; 
  TString outputFilename3 = "Protons.QA.Histograms.root"; 
  TString outputFilename4 = "Protons.Efficiency.root"; 

  printf("****** Connect to PROOF *******\n");
  TProof::Open("alicecaf.cern.ch"); 
  gProof->SetParallel();

  // Enable the Analysis Package
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  gProof->UploadPackage("CORRFW.par");
  gProof->EnablePackage("CORRFW");
  gProof->UploadPackage("PWG2spectra.par");
  gProof->EnablePackage("PWG2spectra");
  
  gProof->Load("AliAnalysisTaskProtonsQA.cxx++");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtonsQA *taskProtonsQA = new AliAnalysisTaskProtonsQA("TaskProtonsQA");
  taskProtonsQA->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
  taskProtonsQA->SetAcceptedVertexDiamond(5.,5.,15.);
  mgr->AddTask(taskProtonsQA);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
							   TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("globalQAList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename1.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("pdgCodeList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename2.Data());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("mcProcessList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename2.Data());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("acceptedCutList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename3.Data());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("acceptedDCAList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename3.Data());
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("efficiencyList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename4.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  if(dataset)
    mgr->StartAnalysis("proof",dataset,stats);
  else {
    // You should get this macro and the txt file from:
    // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
    gROOT->LoadMacro("CreateESDChain.C");
    TChain* chain = 0x0;
    chain = CreateESDChain("ESD82XX_30K.txt",stats);
    chain->SetBranchStatus("*Calo*",0);

    mgr->StartAnalysis("proof",chain);
    //mgr->StartAnalysis("local",chain);
  }
}

//_________________________________________________//
Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    //gROOT->ProcessLine(processline);
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
