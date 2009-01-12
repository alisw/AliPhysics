void runProtonAnalysis(const char* esdAnalysisType = "TPC",) {
  //Macro to run the proton analysis tested for local, proof & GRID.
  //Local: Takes three arguments, the analysis mode, the type of the ESD 
  //       analysis and the path where the tag and ESD or AOD files reside.
  //Interactive: Takes three arguments, the analysis mode, the type of the ESD 
  //             analysis and the name of the collection of tag files.
  //Batch: Takes three arguments, the analysis mode, the type of the ESD 
  //       analysis and the name of the collection file with the event list 
  //       for each file.
  //Proof: Takes four arguments, the analysis mode, the number of events,
  //       the dataset name and the analysis type in case of ESD.
  
  //Analysis mode can be: "MC", "ESD", "AOD"
  //ESD analysis type can be one of the three: "TPC", "Hybrid", "Global"
  TStopwatch timer;
  timer.Start();
  
  //runLocal("ESD",esdAnalysisType,"/home/pchrist/ALICE/Alien/Tutorial/November2007/Tags");
  //runInteractive("ESD",esdAnalysisType,"tag.xml");
  //runBatch("ESD",esdAnalysisType,"wn.xml");  
  runProof("MC",
	   200000,
	   "/COMMON/COMMON/LHC08c11_10TeV_0.5T",
	   esdAnalysisType);
  
  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* mode = "ESD",
	      const char* analysisType = 0x0,
	      const char* path = "/home/pchrist/ALICE/Alien/Tutorial/November2007/Tags") {
  TStopwatch timer;
  timer.Start();

  TString smode = mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
  outputFilename += ".root";

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
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  //____________________________________________________//  

  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD"); 
  tagAnalysis->ChainLocalTags(path);

  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  
  TChain* chain = 0x0;
  chain = tagAnalysis->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
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
  taskProtons->SetTriggerMode(AliAnalysisTaskProtons::kMB2);
  switch(analysisType) {
  case "TPC":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kTPC);
    break;
  case "Hybrid":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kHybrid);
    break;
  case "Global":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kGlobal);
    break;
  default:
    break;
  }
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
                                                            outputFilename.Data());

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
void runInteractive(const char* mode = "ESD",
		    const char* analysisType = 0x0,
		    const char* collectionName = "tag.xml") {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libProofPlayer.so");

  TString smode = mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
  outputFilename += ".root";

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
 
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
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  //____________________________________________________//  
  
  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
 
  //grid tags
  TAlienCollection* coll = TAlienCollection::Open(collectionName);
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  tagAnalysis->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = tagAnalysis->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);
  
  //____________________________________________//
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
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
  taskProtons->SetTriggerMode(AliAnalysisTaskProtons::kMB2);
  switch(analysisType) {
  case "TPC":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kTPC);
    break;
  case "Hybrid":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kHybrid);
    break;
  case "Global":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kGlobal);
    break;
  default:
    break;
  }
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
                                                            outputFilename.Data());
  
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
void runBatch(const char* mode = "ESD",
	      const char* analysisType = 0x0,
	      const char *collectionfile = "wn.xml") {
  TStopwatch timer;
  timer.Start();

  TString smode = mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
  outputFilename += ".root";

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");

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
  setupPar->UploadPackage("CORRFW.par");
  gSystem->EnablePackage("CORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra.so");
  //____________________________________________________//  

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = tagAnalysis->GetChainFromCollection(collectionfile,"esdTree");
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  gROOT->LoadMacro("AliAnalysisTaskProtons.cxx++");
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
  taskProtons->SetTriggerMode(AliAnalysisTaskProtons::kMB2);
  switch(analysisType) {
  case "TPC":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kTPC);
    break;
  case "Hybrid":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kHybrid);
    break;
  case "Global":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kGlobal);
    break;
  default:
    break;
  }
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
                                                            outputFilename.Data());

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
	      const char* dataset = 0x0,
	      const char* analysisType = 0x0) {
  TStopwatch timer;
  timer.Start();
  
  TString smode = mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
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
  
  //____________________________________________//
  gProof->Load("AliAnalysisTaskProtons.cxx++");
  //____________________________________________//

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
  taskProtons->SetTriggerMode(AliAnalysisTaskProtons::kMB2);
  switch(analysisType) {
  case "TPC":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kTPC);
    break;
  case "Hybrid":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kHybrid);
    break;
  case "Global":
    taskProtons->SetAnalysisMode(AliAnalysisTaskProtons::kGlobal);
    break;
  default:
    break;
  }
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
