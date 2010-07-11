void runProtonAnalysis(Bool_t kAnalyzeMC = kTRUE,
		       const char* esdAnalysisType = "Hybrid",
		       const char* pidMode = "Ratio",
		       Bool_t kUseOnlineTrigger = kTRUE,
		       Bool_t kUseOfflineTrigger = kTRUE,
		       Bool_t kRunQA = kFALSE) {
  //Macro to run the proton analysis tested for local, proof & GRID.
  //Local: Takes six arguments, the analysis mode, a boolean to define the ESD
  //       analysis of MC data, the type of the ESD analysis, the PID mode, 
  //       the run number for the offline trigger in case of real data 
  //       analysis and the path where the tag and ESD or AOD files reside.
  //Interactive: Takes six arguments, the analysis mode, a boolean to define 
  //             the ESD analysis of MC data, the type of the ESD analysis, 
  //             the PID mode, the run number for the offline trigger in case 
  //             of real data analysis and the name of the collection of tag 
  //             files.
  //Batch: Takes six arguments, the analysis mode, a boolean to define 
  //       the ESD analysis of MC data, the type of the ESD analysis, 
  //       the PID mode, the run number for the offline trigger in case 
  //       of real data analysis and the name of the collection file with 
  //       the event list for each file.
  //Proof: Takes eight arguments, the analysis mode, a boolean to define 
  //       the ESD analysis of MC data, the type of the ESD analysis, 
  //       the PID mode, the run number for the offline trigger in case 
  //       of real data analysis, the number of events to be analyzed, 
  //       the event number from where we start the analysis and the dataset 
  //========================================================================
  //Analysis mode can be: "MC", "ESD", "AOD"
  //ESD analysis type can be one of the three: "TPC", "Hybrid", "Global"
  //PID mode can be one of the four: "Bayesian" (standard Bayesian approach) 
  //   "Ratio" (ratio of measured over expected/theoretical dE/dx a la STAR) 
  //   "Sigma" (N-sigma area around the fitted dE/dx vs P band)
  TStopwatch timer;
  timer.Start();
  
  /*runLocal("ESD", 
	   kAnalyzeMC,
	   esdAnalysisType,
	   pidMode, kUseOnlineTrigger,kUseOfflineTrigger,
	   "/home/pchrist/ALICE/Baryons/Data/104070");*/
  //runInteractive("ESD", kAnalyzeMC, esdAnalysisType, pidMode, kUseOnlineTrigger, kUseOfflineTrigger, "tag.xml");
  //runBatch("ESD", kAnalyzeMC, esdAnalysisType, pidMode, kUseOnlineTrigger, kUseOfflineTrigger, "wn.xml");  
  runProof("ESD", kAnalyzeMC, esdAnalysisType, pidMode, kUseOnlineTrigger, 
	   kUseOfflineTrigger,
	   500000,0,"/COMMON/COMMON/LHC10a8_run104867_8#esdTree");

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* mode = "ESD",
	      Bool_t kAnalyzeMC = kTRUE,
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0,
	      Bool_t kUseOnlineTrigger = kTRUE,
	      Bool_t kUseOfflineTrigger = kTRUE,
	      const char* path = "/home/pchrist/ALICE/Alien/Tutorial/November2007/Tags") {
  TString smode = mode;
  TString cutFilename = "ListOfCuts."; cutFilename += mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    cutFilename += "."; cutFilename += analysisType;
    outputFilename += "."; outputFilename += analysisType;
  }
  if(pidMode) {
    cutFilename += "."; cutFilename += pidMode;
    outputFilename += "."; outputFilename += pidMode;
  }
 cutFilename += ".root";
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
  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
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
  gROOT->LoadMacro("configProtonAnalysis.C");
  AliProtonAnalysis *analysis = GetProtonAnalysisObject(mode,kAnalyzeMC,
							analysisType,
							pidMode,
							kUseOnlineTrigger,
							kUseOfflineTrigger,
							kRunQA);
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  
  //____________________________________________//
  //Create the proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetAnalysisObject(analysis);
  mgr->AddTask(taskProtons);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("outputQAList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  /*AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cutCanvas",
                                                            TCanvas::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());*/

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  mgr->ConnectOutput(taskProtons,1,coutput2);
  //mgr->ConnectOutput(taskProtons,2,coutput3);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runInteractive(const char* mode = "ESD",
		    Bool_t kAnalyzeMC = kTRUE,
		    const char* analysisType = 0x0,
		    const char* pidMode = 0x0,
		    Bool_t kUseOnlineTrigger = kTRUE,
		    Bool_t kUseOfflineTrigger = kTRUE,
		    const char* collectionName = "tag.xml") {
  gSystem->Load("libProofPlayer.so");
  
  TString smode = mode;
  TString cutFilename = "ListOfCuts."; cutFilename += mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    cutFilename += "."; cutFilename += analysisType;
    outputFilename += "."; outputFilename += analysisType;
  }
  if(pidMode) {
    cutFilename += "."; cutFilename += pidMode;
    outputFilename += "."; outputFilename += pidMode;
  }
  cutFilename += ".root";
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
  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
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
  gROOT->LoadMacro("configProtonAnalysis.C");
  AliProtonAnalysis *analysis = GetProtonAnalysisObject(mode,kAnalyzeMC,
							analysisType,
							pidMode,
							kUseOnlineTrigger,
							kUseOfflineTrigger,
							kRunQA);
  //runNumberForOfflineTtrigger);
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }

  //____________________________________________//
  //Create the proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetAnalysisObject(analysis);
  mgr->AddTask(taskProtons);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("outputQAList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  /*AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cutCanvas",
                                                            TCanvas::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());*/

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  mgr->ConnectOutput(taskProtons,1,coutput2);
  //mgr->ConnectOutput(taskProtons,2,coutput3);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runBatch(const char* mode = "ESD",
	      Bool_t kAnalyzeMC = kTRUE,
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0,
	      Bool_t kUseOnlineTrigger = kTRUE,
	      Bool_t kUseOfflineTrigger = kTRUE,
	      const char *collectionfile = "wn.xml") {
  TString smode = mode;
  TString cutFilename = "ListOfCuts."; cutFilename += mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    cutFilename += "."; cutFilename += analysisType;
    outputFilename += "."; outputFilename += analysisType;
  }
  if(pidMode) {
    cutFilename += "."; cutFilename += pidMode;
    outputFilename += "."; outputFilename += pidMode;
  }
  cutFilename += ".root";
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
  setupPar("CORRFW");
  gSystem->Load("libCORRFW.so");
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
  gROOT->LoadMacro("configProtonAnalysis.C");
  AliProtonAnalysis *analysis = GetProtonAnalysisObject(mode,kAnalyzeMC,
							analysisType,
							pidMode,
							kUseOnlineTrigger,
							kUseOfflineTrigger,
							kRunQA);
  //runNumberForOfflineTtrigger);
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  
  //____________________________________________//
  //Create the proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetAnalysisObject(analysis);
  mgr->AddTask(taskProtons);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("outputQAList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  /*AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cutCanvas",
                                                            TCanvas::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());*/
  
  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  mgr->ConnectOutput(taskProtons,1,coutput2);
  //mgr->ConnectOutput(taskProtons,2,coutput3);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid",chain);
}

//_________________________________________________//
void runProof(const char* mode = "ESD",
	      Bool_t kAnalyzeMC = kTRUE,
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0,
	      Bool_t kUseOnlineTrigger = kTRUE,
	      Bool_t kUseOfflineTrigger = kTRUE,
	      Int_t stats = 0, Int_t startingPoint = 0,
	      const char* dataset = 0x0) {  
  TString smode = mode;
  TString cutFilename = "ListOfCuts."; cutFilename += mode;
  TString outputFilename = "Protons."; outputFilename += mode;
  if(analysisType) {
    cutFilename += "."; cutFilename += analysisType;
    outputFilename += "."; outputFilename += analysisType;
  }
  if(pidMode) {
    cutFilename += "."; cutFilename += pidMode;
    outputFilename += "."; outputFilename += pidMode;
  }
  cutFilename += ".root";
  outputFilename += ".root";

  gEnv->SetValue("XSec.GSI.DelegProxy","2");
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
  gROOT->LoadMacro("configProtonAnalysis.C");
  AliProtonAnalysis *analysis = GetProtonAnalysisObject(mode,kAnalyzeMC,
							analysisType,
							pidMode,
							kUseOnlineTrigger,
							kUseOfflineTrigger,
							kRunQA);
  //____________________________________________//

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  if(smode == "MC") {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  //____________________________________________//
  //Create the proton task
  AliAnalysisTaskProtons *taskProtons = new AliAnalysisTaskProtons("TaskProtons");
  taskProtons->SetAnalysisObject(analysis);
  mgr->AddTask(taskProtons);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("outputQAList",
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());
  /*AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cutCanvas",
                                                            TCanvas::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            outputFilename.Data());*/

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  mgr->ConnectOutput(taskProtons,1,coutput2);
  //mgr->ConnectOutput(taskProtons,3,coutput3);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  if(dataset)
    mgr->StartAnalysis("proof",dataset,stats,startingPoint);
  else {
    // You should get this macro and the txt file from:
    // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
    gROOT->LoadMacro("CreateESDChain.C");
    TChain* chain = 0x0;
    chain = CreateESDChain("ESD82XX_30K.txt",stats);
    chain->SetBranchStatus("*Calo*",0);

    mgr->StartAnalysis("proof",chain);
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
