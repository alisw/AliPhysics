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
    kRunQA,
    "/home/pchrist/ALICE/Baryons/Data/104070");*/
  //runInteractive("ESD", kAnalyzeMC, esdAnalysisType, pidMode, kUseOnlineTrigger, kUseOfflineTrigger, kRunQA, "tag.xml");
  //runBatch("ESD", kAnalyzeMC, esdAnalysisType, pidMode, kUseOnlineTrigger, kUseOfflineTrigger, kRunQA, "wn.xml");  
  runLocal("ESD", kAnalyzeMC, esdAnalysisType, pidMode,
           kUseOnlineTrigger, kUseOfflineTrigger, kRunQA, ".");

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
	      Bool_t kRunQA = kFALSE,
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
  gSystem->Load("libSTEERBase");
  setupPar("ESD");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  setupPar("AOD");
  gSystem->Load("libAOD");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice");
  setupPar("CORRFW");
  gSystem->Load("libCORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra");
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
		    Bool_t kRunQA = kFALSE,
		    const char* collectionName = "tag.xml") {
  
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
  gSystem->Load("libSTEERBase");
  setupPar("ESD");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  setupPar("AOD");
  gSystem->Load("libAOD");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice");
  setupPar("CORRFW");
  gSystem->Load("libCORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra");
  //____________________________________________________//  
  
  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
 
  //grid tags
  TGridCollection* coll = gGrid->OpenCollection(collectionName);
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
	      Bool_t kRunQA = kFALSE,
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

  //____________________________________________________//
  //_____________Setting up the par files_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase");
  setupPar("ESD");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  setupPar("AOD");
  gSystem->Load("libAOD");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice");
  setupPar("CORRFW");
  gSystem->Load("libCORRFW");
  setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra");
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
