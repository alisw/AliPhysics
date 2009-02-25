void runProtonAnalysisQA(const char *analysisType = "TPC") {
  //Macro to run the proton QA analysis tested for local, proof & GRID.
  //Proof: Takes three arguments, the number of events, the dataset name 
  //       and the analysis type
  //Interactive: Takes two arguments, the name of the collection file with 
  //             the event list for each file and the analysis type
  //Batch: Takes two arguments, the name of the collection file with 
  //             the ESDs and the analysis type
  //Analysis type can be one of the three: "TPC", "Hybrid", "Global"
  TStopwatch timer;
  timer.Start();
  
  //runProof(200000,"/COMMON/COMMON/LHC08c11_10TeV_0.5T",analysisType);
  //runInteractive("wn.xml",analysisType);
  //runBatch("wn.xml",analysisType);
  runLocal(analysisType);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* analysisType) {
  TString outputFilename1 = "Protons.QA."; outputFilename1 += analysisType;
  outputFilename1 += ".root"; //main QA file
  TString outputFilename2 = "Protons.MC.QA."; outputFilename2 += analysisType;
  outputFilename2 += ".root"; //MC process QA
  TString outputFilename3 = "Protons.QA.Histograms."; 
  outputFilename3 += analysisType;
  outputFilename3 += ".root"; //Accepted cut distributions
  TString outputFilename4 = "Protons.Efficiency."; 
  outputFilename4 += analysisType;
  outputFilename4 += ".root"; //Reco and PID efficiency
  TString outputFilename5 = "Vertex.QA.root"; //vertex QA

  gSystem->Load("libProofPlayer.so");

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
  //Create the chain from the xml collection of esd files
  TChain *chain = new TChain("esdTree");  
  chain->Add("AliESDs.root");
  
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
  switch(analysisType) {
  case "TPC":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
    break;
  case "Hybrid":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
    break;
  case "Global":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kGlobal);
    break;
  default:
    break;
  }
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
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("vertexList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename5.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  mgr->ConnectOutput(taskProtonsQA,6,coutput7);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runBatch(const char *collectionfile,
	      const char *analysisType) {
  TString outputFilename1 = "Protons.QA."; outputFilename1 += analysisType;
  outputFilename1 += ".root"; //main QA file
  TString outputFilename2 = "Protons.MC.QA."; outputFilename2 += analysisType;
  outputFilename2 += ".root"; //MC process QA
  TString outputFilename3 = "Protons.QA.Histograms."; 
  outputFilename3 += analysisType;
  outputFilename3 += ".root"; //Accepted cut distributions
  TString outputFilename4 = "Protons.Efficiency."; 
  outputFilename4 += analysisType;
  outputFilename4 += ".root"; //Reco and PID efficiency
  TString outputFilename5 = "Vertex.QA.root"; //vertex QA

  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");

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
  //Create the chain from the xml collection of esd files
  TAlienCollection* coll = TAlienCollection::Open(collectionfile);
  TGridResult* gridResult = coll->GetGridResult("",0,0);
  Int_t nEntries = gridResult->GetEntries();
  TChain *chain = new TChain("esdTree");  
  TString alienUrl = 0;
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = gridResult->GetKey(i,"turl");
    chain->Add(alienUrl);
  }
  
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
  switch(analysisType) {
  case "TPC":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
    break;
  case "Hybrid":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
    break;
  case "Global":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kGlobal);
    break;
  default:
    break;
  }
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
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("vertexList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename5.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  mgr->ConnectOutput(taskProtonsQA,6,coutput7);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runInteractive(const char *collectionfile,
		    const char *analysisType) {
  TString outputFilename1 = "Protons.QA."; outputFilename1 += analysisType;
  outputFilename1 += ".root"; //main QA file
  TString outputFilename2 = "Protons.MC.QA."; outputFilename2 += analysisType;
  outputFilename2 += ".root"; //MC process QA
  TString outputFilename3 = "Protons.QA.Histograms."; 
  outputFilename3 += analysisType;
  outputFilename3 += ".root"; //Accepted cut distributions
  TString outputFilename4 = "Protons.Efficiency."; 
  outputFilename4 += analysisType;
  outputFilename4 += ".root"; //Reco and PID efficiency
  TString outputFilename5 = "Vertex.QA.root"; //vertex QA

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
  switch(analysisType) {
  case "TPC":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
    break;
  case "Hybrid":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
    break;
  case "Global":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kGlobal);
    break;
  default:
    break;
  }
  taskProtonsQA->SetAcceptedVertexDiamond(5.,5.,15.);
  mgr->AddTask(taskProtonsQA);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
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
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("vertexList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename5.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  mgr->ConnectOutput(taskProtonsQA,6,coutput7);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}
 
//_________________________________________________//
void runProof(Int_t stats = 0, 
	      const char* dataset = 0x0,
	      const char *analysisType) {
  TString outputFilename1 = "Protons.QA."; outputFilename1 += analysisType;
  outputFilename1 += ".root"; //main QA file
  TString outputFilename2 = "Protons.MC.QA."; outputFilename2 += analysisType;
  outputFilename2 += ".root"; //MC process QA
  TString outputFilename3 = "Protons.QA.Histograms."; 
  outputFilename3 += analysisType;
  outputFilename3 += ".root"; //Accepted cut distributions
  TString outputFilename4 = "Protons.Efficiency."; 
  outputFilename4 += analysisType;
  outputFilename4 += ".root"; //Reco and PID efficiency
  TString outputFilename5 = "Vertex.QA.root"; //vertex QA

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
  switch(analysisType) {
  case "TPC":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
    break;
  case "Hybrid":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kHybrid);
    break;
  case "Global":
    taskProtonsQA->SetAnalysisMode(AliAnalysisTaskProtonsQA::kGlobal);
    break;
  default:
    break;
  }
  taskProtonsQA->SetAcceptedVertexDiamond(5.,5.,15.);
  mgr->AddTask(taskProtonsQA);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
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
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("vertexList", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename5.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtonsQA,0,cinput1);
  mgr->ConnectOutput(taskProtonsQA,0,coutput1);
  mgr->ConnectOutput(taskProtonsQA,1,coutput2);
  mgr->ConnectOutput(taskProtonsQA,2,coutput3);
  mgr->ConnectOutput(taskProtonsQA,3,coutput4);
  mgr->ConnectOutput(taskProtonsQA,4,coutput5);
  mgr->ConnectOutput(taskProtonsQA,5,coutput6);
  mgr->ConnectOutput(taskProtonsQA,6,coutput7);
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
