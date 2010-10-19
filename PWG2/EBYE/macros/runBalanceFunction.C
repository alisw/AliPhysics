enum anaModes {kLocal,kLocalPAR,kPROOF,kInteractive,kGRID};
//kLocal: Analyze locally files in your computer using aliroot
//kLocalPAR: Analyze locally files in your computer using root + PAR files
//kPROOF: Analyze CAF files with PROOF

void runBalanceFunction(Int_t analysisMode = kInteractive,
			Bool_t kMCAnalysis = kFALSE,
			const char* dataMode = "ESD") {
  TStopwatch timer;
  timer.Start();

  //Load the libraries
  LoadLibraries(analysisMode);

  //Select the running mode and create the chain
  if (analysisMode == kGRID) {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler();  
    if (!alienHandler) return;
  }
  if (analysisMode==kLocal || analysisMode == kLocalPAR) {
    TChain *chain = new TChain("esdTree");
    chain->Add("../Set1/AliESDs.root");
    chain->Add("../Set2/AliESDs.root");
  }
  if(analysisMode == kInteractive) {
    TGrid::Connect("alien://");
    TChain *chain = new TChain("esdTree");
    TString alienUrl;

    TAlienCollection *collection = TAlienCollection::Open("wn.xml");
    TGridResult *gResult = collection->GetGridResult("",0,0);
    Int_t nEntries = gResult->GetEntries();
    for(Int_t i = 0; i < nEntries; i++) {
      alienUrl = gResult->GetKey(i,"turl");
      chain->Add(alienUrl.Data());
    }
  }
  //___________________________________________________//
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  if (analysisMode == kGRID) { 
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  //Configure the BF object
  gROOT->LoadMacro("AddTaskBalanceFunction.C");
  AliAnalysisTaskBF *taskBF = AddTaskBalanceFunction();

  // Task to check the offline trigger
  if (analysisMode == kLocal || analysisMode == kGRID || analysisMode == kInteractive) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C"); }
  else if (analysisMode == kPROOF || analysisMode == kLocalPAR) {
    gROOT->LoadMacro("AddTaskPhysicsSelection.C"); }
  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();
  if (kMCAnalysis) {physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();}

  // Enable debug printouts
  mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  
  if (analysisMode == kLocal || analysisMode == kLocalPAR || analysisMode == kInteractive) {
    mgr->StartAnalysis("local",chain);
  }
  else if (analysisMode == kPROOF) {
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  }
  else if (analysisMode == kGRID) { 
    mgr->StartAnalysis("grid");
  }

  timer.Stop();
  timer.Print();
}

void SetupPar(char* pararchivename) {
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
  
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 
  if ( gSystem->AccessPathName(pararchivename) ) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
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
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

//__________________________________________________________//
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

//__________________________________________________________//
void LoadLibraries(const anaModes mode) {  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==kLocal || mode==kGRID || mode==kInteractive) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    if (mode==kLocal) {
      Printf("Local: loading the libPWG2ebye.so");
      gSystem->Load("libPWG2ebye.so");
    }
    if (mode==kGRID || mode==kInteractive) {
      //setupPar("PWG2ebye");
      //Printf("GRID: loading the libPWG2ebye.so");
      gSystem->Load("libPWG2ebye.so");
    }
  }//local or GRID
  else if (mode == kLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    gSystem->Load("libSTEERBase.so");
    //setupPar("ESD");
    gSystem->Load("libVMC.so");
    gSystem->Load("libESD.so");
    //setupPar("AOD");
    gSystem->Load("libAOD.so");
    //setupPar("ANALYSIS");
    gSystem->Load("libANALYSIS.so");
    //setupPar("ANALYSISalice");
    gSystem->Load("libANALYSISalice.so");
    Int_t setuparflag = setupPar("PWG2ebye");
    Printf("localPar: loading the libPWG2ebye.so (%d)",setuparflag);
    gSystem->Load("libPWG2ebye.so");
  }//local with par files
  
  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==kPROOF) {
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof::Open("alice-caf.cern.ch");
 
    // Upload the Packages
    gProof->UploadPackage("STEERBase.par");
    gProof->UploadPackage("ESD.par");    
    gProof->UploadPackage("AOD.par");       
    gProof->UploadPackage("ANALYSIS.par"); 
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->UploadPackage("CORRFW.par");
    gProof->UploadPackage("PWG2ebye.par");

    // Enable the Packages     
    gProof->EnablePackage("STEERBase");
    gProof->EnablePackage("ESD");
    gProof->EnablePackage("AOD");
    gProof->EnablePackage("ANALYSIS");
    gProof->EnablePackage("ANALYSISalice");
    gProof->EnablePackage("PWG2ebye");

    // Show enables Packages
    gProof->ShowEnabledPackages();
  }  
  
}
