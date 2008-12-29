void runProtonAnalysis() {
  TStopwatch timer;
  timer.Start();
  
  //runLocal("ESD");
  //runInteractive("ESD");
  //runBatch("ESD");
  
  runProof("ESD",200000,"/COMMON/COMMON/LHC08c11_10TeV_0.5T"); //use data sets
  //runProof("ESD",200); //use ascii files
  
  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* mode = "ESD") {
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

  //_________________________________________________________//
  //___________Setting up ANALYSISalice.par__________________//
  //_________________________________________________________//
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  //__________________________________________________//
  //___________Setting up CORRFW.par__________________//
  //__________________________________________________//
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");

  //____________________________________________________________//
  //_____________Setting up PWG2spectra.par_____________________//
  //____________________________________________________________//
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");

  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD"); 
  tagAnalysis->ChainLocalTags("/home/pchrist/ALICE/Alien/Tutorial/November2007/Tags");

  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  
  TChain* chain = 0x0;
  chain = tagAnalysis->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }

  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetType(mode);
  taskProtons->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtons->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
  taskProtons->SetAcceptedVertexDiamond(5.,5.,15.);
  //Momentum dependent priors
  /*TFile *f = TFile::Open("PriorProb/PriorProbabilities.root ");
  TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
  TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
  TF1 *fitPions = (TF1 *)f->Get("fitPions");
  TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
  TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
  taskProtons->SetPriorProbabilityFunctions(fitElectrons,
					    fitMuons,
					    fitPions,
					    fitKaons,
					    fitProtons);*/
  mgr->AddTask(taskProtons);

  // Create containers for input/output                                                                              
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
                                                           TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList1",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputCont
                                                            "Protons.ESD.root");

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runInteractive(const char* mode = "ESD") {
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

  //_________________________________________________________//
  //___________Setting up ANALYSISalice.par__________________//
  //_________________________________________________________//
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  //__________________________________________________//
  //___________Setting up CORRFW.par__________________//
  //__________________________________________________//
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");

  //____________________________________________________________//
  //_____________Setting up PWG2spectra.par_____________________//
  //____________________________________________________________//
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
  
  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
 
  //grid tags
  TAlienCollection* coll = TAlienCollection::Open("tag.xml");
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  tagAnalysis->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = tagAnalysis->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }

  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetType(mode);
  taskProtons->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtons->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
  taskProtons->SetAcceptedVertexDiamond(5.,5.,15.);
  //Momentum dependent priors
  /*TFile *f = TFile::Open("PriorProb/PriorProbabilities.root ");
  TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
  TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
  TF1 *fitPions = (TF1 *)f->Get("fitPions");
  TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
  TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
  taskProtons->SetPriorProbabilityFunctions(fitElectrons,
					    fitMuons,
					    fitPions,
					    fitKaons,
					    fitProtons);*/
  mgr->AddTask(taskProtons);

  // Create containers for input/output                                                                               
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
                                                           TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList1",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputCont
                                                            "Protons.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runBatch(const char* mode = "ESD") {
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

  //_________________________________________________________//
  //___________Setting up ANALYSISalice.par__________________//
  //_________________________________________________________//
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");

  //__________________________________________________//
  //___________Setting up CORRFW.par__________________//
  //__________________________________________________//
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");

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
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = tagAnalysis->GetChainFromCollection(collectionfile,"esdTree");
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetType(mode);
  taskProtons->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtons->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
  taskProtons->SetAcceptedVertexDiamond(5.,5.,15.);
  //Momentum dependent priors
  /*TFile *f = TFile::Open("PriorProb/PriorProbabilities.root ");
  TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
  TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
  TF1 *fitPions = (TF1 *)f->Get("fitPions");
  TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
  TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
  taskProtons->SetPriorProbabilityFunctions(fitElectrons,
					    fitMuons,
					    fitPions,
					    fitKaons,
					    fitProtons);*/
  mgr->AddTask(taskProtons);

  // Create containers for input/output                                                                               
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
                                                           TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList1",
                                                            TList::Class(),AliAnalysisManager::kOutputCont
                                                            "Protons.ESD.root");

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid",chain);

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runProof(const char* mode = "ESD", 
	      Int_t stats = 0, 
	      const char* dataset = 0x0) {
  TStopwatch timer;
  timer.Start();
  
  TString smode = mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  outputFilename += ".root";

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
  
  gProof->Load("AliAnalysisTaskProtons.cxx++");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetType(mode);
  taskProtons->SetTriggerMode(AliAnalysisTaskProtonsQA::kMB2);
  taskProtons->SetAnalysisMode(AliAnalysisTaskProtonsQA::kTPC);
  taskProtons->SetAcceptedVertexDiamond(5.,5.,15.);
  //Momentum dependent priors
  /*TFile *f = TFile::Open("PriorProb/PriorProbabilities.root ");
  TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
  TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
  TF1 *fitPions = (TF1 *)f->Get("fitPions");
  TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
  TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
  taskProtons->SetPriorProbabilityFunctions(fitElectrons,
					    fitMuons,
					    fitPions,
					    fitKaons,
					    fitProtons);*/
  mgr->AddTask(taskProtons);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("dataChain",
							   TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList1", 
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFilename.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
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
  }

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
