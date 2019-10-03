void runProtonsCorrectionAnalysis(const char* analysisType  = "Hybrid",
				  const char* pidMode = "Ratio",
				  Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, 
				  Bool_t fIsOn_AliProtonFeedDownAnalysis=kFALSE,
				  Bool_t fIsOn_AliProtonSpectraCorrection=kFALSE) 
{
  //Macro to run the proton feed-down analysis tested for local, proof & GRID.
  //Local: Takes four arguments, the analysis mode, the type of the ESD 
  //       analysis, the PID mode and the path where the tag and ESD or 
  //       AOD files reside.
  //Interactive: Takes four arguments, the analysis mode, the type of the ESD 
  //             analysis, the PID mode and the name of the collection of tag 
  //             files.
  //Batch: Takes four arguments, the analysis mode, the type of the ESD 
  //       analysis, the PID mode and the name of the collection file with 
  //       the event list for each file.
  //Proof: Takes five arguments, the analysis level, the analysis mode in 
  //       case of ESD, the PID mode, the number of events and the dataset 
  //       name and .  
  //Analysis mode can be: "MC", "ESD", "AOD"
  //ESD analysis type can be one of the three: "TPC", "Hybrid", "Global"
  //PID mode can be one of the four: "Bayesian" (standard Bayesian approach) 
  //   "Ratio" (ratio of measured over expected/theoretical dE/dx a la STAR) 
  //   "Sigma1" (N-sigma area around the fitted dE/dx vs P band)
  //   "Sigma2" (same as previous but taking into account the No of TPC points)
  TStopwatch timer;
  timer.Start();
  
  //runLocal("ESD",analysisType,pidMode,"/home/pchrist/ALICE/Baryons/Analysis/Protons/Local/data",fIsOn_AliProtonAbsorptionCorrection,fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection);
  //runInteractive("ESD",analysisType,pidMode,"/home/marek/Analysis/global_xml/tagtest81627.xml",fIsOn_AliProtonAbsorptionCorrection,fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection);
  //runBatch("ESD",analysisType,pidMode,"wn.xml",fIsOn_AliProtonAbsorptionCorrection,fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection);  
  runProof("ESD",analysisType,pidMode,500000,"/COMMON/COMMON/LHC10a8_run104867_8#esdTree",fIsOn_AliProtonAbsorptionCorrection,fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection);
  
  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* mode = "ESD",
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0, 
	      const char* path = 0x0,
	      Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, 
	      Bool_t fIsOn_AliProtonFeedDownAnalysis=kTRUE,
	      Bool_t fIsOn_AliProtonSpectraCorrection=kTRUE) {
  TString smode = mode;
  TString outputFilename = "ProtonCorrection."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
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
  //chain->SetBranchStatus("*Calo*",0);
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManagerProtonCorrection");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  
  //____________________________________________//
  /*AliProtonCorrectionAnalysisTask *taskProtons = new AliProtonCorrectionAnalysisTask("TaskProtonsProtonCorrection");
    
  if(fIsOn_AliProtonAbsorptionCorrection||fIsOn_AliProtonFeedDownAnalysis||fIsOn_AliProtonSpectraCorrection)
  {
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(mode,analysisType,pidMode);
  taskProtons->SetBaseAnalysis(baseAnalysis);
  }	
  else
  return;
  if(fIsOn_AliProtonAbsorptionCorrection)
  {
  AliProtonAbsorptionCorrection * absorptioncorrection=new AliProtonAbsorptionCorrection();
  taskProtons->SetAnalysisObjectAbsorptionCorrection(absorptioncorrection);
  }
  if(fIsOn_AliProtonFeedDownAnalysis)
  {
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
  AliProtonFeedDownAnalysis *analysisFeedDown = GetProtonFeedDownAnalysisObject();
  taskProtons->SetAnalysisObjectFeedDown(analysisFeedDown);
  }	
  if(fIsOn_AliProtonSpectraCorrection)
  {
  AliProtonSpectraCorrection* spectracorrection=new AliProtonSpectraCorrection();
  taskProtons->SetAnalysisObjectSpectraCorrection(spectracorrection);
  }*/
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonCorrectionAnalysisTask.C"); 
  AliProtonCorrectionAnalysisTask* taskProtons=GetAliProtonCorrectionAnalysisTask(mode,analysisType,pidMode ,fIsOn_AliProtonAbsorptionCorrection, fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection); 
  mgr->AddTask(taskProtons);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",TList::Class(),AliAnalysisManager::kOutputContainer,outputFilename.Data());
  
  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runInteractive(const char* mode = "ESD",
		    const char* analysisType = 0x0,
		    const char* pidMode = 0x0,
		    const char* collectionName = "tag.xml",
		    Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, 
		    Bool_t fIsOn_AliProtonFeedDownAnalysis=kTRUE,
		    Bool_t fIsOn_AliProtonSpectraCorrection=kTRUE) {
  gSystem->Load("libProofPlayer");
  
  TString smode = mode;
  TString outputFilename = "ProtonCorrection."; outputFilename += mode;
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
  //chain->SetBranchStatus("*Calo*",0);
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManagerProtonCorrection");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  
  //____________________________________________//
  /*AliProtonCorrectionAnalysisTask *taskProtons = new AliProtonCorrectionAnalysisTask("TaskProtonsProtonCorrection");
    if(fIsOn_AliProtonAbsorptionCorrection||fIsOn_AliProtonFeedDownAnalysis||fIsOn_AliProtonSpectraCorrection)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
    AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(mode,analysisType,pidMode);
    taskProtons->SetBaseAnalysis(baseAnalysis);
    }	
    else
    return;
    if(fIsOn_AliProtonAbsorptionCorrection)
    {
    AliProtonAbsorptionCorrection* absorptioncorrection=new AliProtonAbsorptionCorrection();
    taskProtons->SetAnalysisObjectAbsorptionCorrection(absorptioncorrection);
    }
    if(fIsOn_AliProtonFeedDownAnalysis)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
    AliProtonFeedDownAnalysis* analysisFeedDown = GetProtonFeedDownAnalysisObject();
    taskProtons->SetAnalysisObjectFeedDown(analysisFeedDown);
    }	
    if(fIsOn_AliProtonSpectraCorrection)
    {
    AliProtonSpectraCorrection* spectracorrection=new AliProtonSpectraCorrection();
    taskProtons->SetAnalysisObjectSpectraCorrection(spectracorrection);
    }*/
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonCorrectionAnalysisTask.C"); 
  AliProtonCorrectionAnalysisTask* taskProtons=GetAliProtonCorrectionAnalysisTask(mode,analysisType,pidMode ,fIsOn_AliProtonAbsorptionCorrection, fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection); 
  mgr->AddTask(taskProtons);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",TList::Class(),AliAnalysisManager::kOutputContainer,outputFilename.Data());
  
  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput1);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runBatch(const char* mode = "ESD",
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0,
	      const char *collectionfile = "wn.xml",
	      Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, 
	      Bool_t fIsOn_AliProtonFeedDownAnalysis=kTRUE,
	      Bool_t fIsOn_AliProtonSpectraCorrection=kTRUE) {
  TString smode = mode;
  TString outputFilename = "ProtonCorrection."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
  outputFilename += ".root";
  
  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer");
  
  //____________________________________________________//
  //_____________Setting up the par files_______________//
  //____________________________________________________//
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW") ;
  
  //setupPar("PWG2spectra");
  gSystem->Load("libPWG2spectra");
  //____________________________________________________//  
  
  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD");
  TChain *chain = 0x0;
  chain = tagAnalysis->GetChainFromCollection(collectionfile,"esdTree");
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManagerProtonCorrection");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  
  //____________________________________________//
  /*AliProtonCorrectionAnalysisTask *taskProtons = new AliProtonCorrectionAnalysisTask("TaskProtonsProtonCorrection");
    if(fIsOn_AliProtonAbsorptionCorrection||fIsOn_AliProtonFeedDownAnalysis||fIsOn_AliProtonSpectraCorrection)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
    AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(mode,analysisType ,pidMode);
    taskProtons->SetBaseAnalysis(baseAnalysis);
    }	
    else
    return;
    if(fIsOn_AliProtonAbsorptionCorrection)
    {
    AliProtonAbsorptionCorrection * absorptioncorrection=new AliProtonAbsorptionCorrection();
    taskProtons->SetAnalysisObjectAbsorptionCorrection(absorptioncorrection);
    }
    if(fIsOn_AliProtonFeedDownAnalysis)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
    AliProtonFeedDownAnalysis *analysisFeedDown = GetProtonFeedDownAnalysisObject();
    taskProtons->SetAnalysisObjectFeedDown(analysisFeedDown);
    }	
    if(fIsOn_AliProtonSpectraCorrection)
    {
    AliProtonSpectraCorrection* spectracorrection=new AliProtonSpectraCorrection();
    taskProtons->SetAnalysisObjectSpectraCorrection(spectracorrection);
    }*/
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonCorrectionAnalysisTask.C"); 
  AliProtonCorrectionAnalysisTask* taskProtons=GetAliProtonCorrectionAnalysisTask(mode,analysisType,pidMode ,fIsOn_AliProtonAbsorptionCorrection, fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection); 
  mgr->AddTask(taskProtons);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",TList::Class(),AliAnalysisManager::kOutputContainer,outputFilename.Data());

  //____________________________________________//
  mgr->ConnectInput(taskProtons,0,cinput);
  mgr->ConnectOutput(taskProtons,0,coutput1);
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runProof(const char* mode = "ESD",
	      const char* analysisType = 0x0,
	      const char* pidMode = 0x0,
	      Int_t stats = 0, 
	      const char* dataset = 0x0,
	      Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, 
	      Bool_t fIsOn_AliProtonFeedDownAnalysis=kTRUE,
	      Bool_t fIsOn_AliProtonSpectraCorrection=kTRUE) {  
  TString smode = mode;
  TString outputFilename = "ProtonCorrection."; outputFilename += mode;
  if(analysisType) {
    outputFilename += "."; outputFilename += analysisType;
  }
  outputFilename += ".root";
  
  printf("****** Connect to PROOF *******\n");
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("alicecaf.cern.ch"); 
  //gProof->SetParallel();
  
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
  
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("protonAnalysisManagerProtonCorrection");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  //____________________________________________//
  //Create the proton task
  /*AliProtonCorrectionAnalysisTask *taskProtons = new AliProtonCorrectionAnalysisTask("TaskProtonsProtonCorrection");
    if(fIsOn_AliProtonAbsorptionCorrection||fIsOn_AliProtonFeedDownAnalysis||fIsOn_AliProtonSpectraCorrection)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
    AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(mode,analysisType ,pidMode);
    taskProtons->SetBaseAnalysis(baseAnalysis);
    }	
    else
    return
    if(fIsOn_AliProtonAbsorptionCorrection)
    {
    AliProtonAbsorptionCorrection * absorptioncorrection=new AliProtonAbsorptionCorrection();
    taskProtons->SetAnalysisObjectAbsorptionCorrection(absorptioncorrection);
    }
    if(fIsOn_AliProtonFeedDownAnalysis)
    {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
    AliProtonFeedDownAnalysis *analysisFeedDown = GetProtonFeedDownAnalysisObject();
    taskProtons->SetAnalysisObjectFeedDown(analysisFeedDown);
    }	
    if(fIsOn_AliProtonSpectraCorrection)
    {
    AliProtonSpectraCorrection* spectracorrection=new AliProtonSpectraCorrection();
    taskProtons->SetAnalysisObjectSpectraCorrection(spectracorrection);
    }*/
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonCorrectionAnalysisTask.C"); 
  AliProtonCorrectionAnalysisTask* taskProtons=GetAliProtonCorrectionAnalysisTask(mode,analysisType,pidMode ,fIsOn_AliProtonAbsorptionCorrection, fIsOn_AliProtonFeedDownAnalysis, fIsOn_AliProtonSpectraCorrection); 
  mgr->AddTask(taskProtons);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputList",TList::Class(),AliAnalysisManager::kOutputContainer,outputFilename.Data());
  
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
