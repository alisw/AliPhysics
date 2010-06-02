enum anaModes {mLocal,mLocalPAR,mPROOF,mGRID};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF

// RUN SETTINGS

// Flow analysis method can be:(set to kTRUE or kFALSE)
Bool_t MCEP     = kFALSE;  // correlation with Monte Carlo reaction plane
Bool_t SP       = kTRUE;  // scalar product method (similar to eventplane method)
Bool_t GFC      = kTRUE;  // cumulants based on generating function
Bool_t QC       = kTRUE;  // cumulants using Q vectors
Bool_t FQD      = kTRUE;  // fit of the distribution of the Q vector (only integrated v)
Bool_t LYZ1SUM  = kTRUE;  // Lee Yang Zeroes using sum generating function (integrated v)
Bool_t LYZ1PROD = kTRUE;  // Lee Yang Zeroes using product generating function (integrated v)
Bool_t LYZ2SUM  = kFALSE; // Lee Yang Zeroes using sum generating function (second pass differential v)
Bool_t LYZ2PROD = kFALSE; // Lee Yang Zeroes using product generating function (second pass differential v)
Bool_t LYZEP    = kFALSE; // Lee Yang Zeroes Event plane using sum generating function (gives eventplane + weight)
Bool_t MH       = kTRUE;  // azimuthal correlators in mixed harmonics  
Bool_t NL       = kTRUE;  // nested loops (for instance distribution of phi1-phi2 for all distinct pairs)

Bool_t METHODS[] = {SP,LYZ1SUM,LYZ1PROD,LYZ2SUM,LYZ2PROD,LYZEP,GFC,QC,FQD,MCEP,MH,NL};

// Analysis type can be ESD, AOD, MC, ESDMCkineESD, ESDMCkineMC
const TString type = "ESD";

// Boolean to fill/not fill the QA histograms
Bool_t QA = kTRUE;   

// Boolean to use/not use weights for the Q vector
Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)


//void runFlowTask(Int_t mode=mLocal, Int_t nRuns = 1, 
//Bool_t DATA = kTRUE, const Char_t* dataDir="/data/alice2/kolk/PP/data/LHC09d/104892/test", Int_t offset = 0)
                //Bool_t DATA = kFALSE, const Char_t* dataDir="/data/alice2/kolk/PP/LHC09d10/104873", Int_t offset = 0)

//void runFlowTask(Int_t mode = mPROOF, Int_t nRuns = 50000000, 
		 //Bool_t DATA = kFALSE, const Char_t* dataDir="/PWG2/akisiel/Therminator_midcentral_ESD", Int_t offset=0)
		 //Bool_t DATA = kFALSE, const Char_t* dataDir="/PWG2/akisiel/LHC10d6_0.9TeV_EPOS_12502X", Int_t offset=0)
		 //Bool_t DATA = kTRUE, const Char_t* dataDir="/alice/data/LHC10b_000115322_p1", Int_t offset=0) // data 7 TeV
		 //Bool_t DATA = kFALSE, const Char_t* dataDir="/alice/sim/LHC10a18_140012", Int_t offset=0) //perugia0 7 TeV
		 //Bool_t DATA = kFALSE, const Char_t* dataDir="/alice/sim/LHC10a20_140514", Int_t offset=0) //phojet 7 TeV

void runFlowTask(Int_t mode = mGRID, Bool_t DATA = kTRUE)
{
  TStopwatch timer;
  timer.Start();
  
  CrossCheckUserSettings(DATA);

  LoadLibraries(mode);

  if (mode == mGRID) {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler();  
    if (!alienHandler) return;
  }
  
  if (mode==mLocal || mode == mLocalPAR) {
    if (type!="AOD") { TChain* chain = CreateESDChain(dataDir, nRuns, offset);}
    else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);}
  }
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
 
  if (mode == mGRID) { 
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
  }

  if (type == "ESD"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
    if (MCEP) { 
      AliMCEventHandler *mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc); 
    }
  }
  
  if (type == "AOD"){
    AliVEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH); 
    if (MCEP) { 
      AliMCEventHandler *mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);
    } 
  }
  
  if (type == "MC" || type == "ESDMCkineESD" || type == "ESDMCkineMC"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
    
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc); 
  }
    
  //____________________________________________//
  // Load the analysis task
  gROOT->LoadMacro("AddTaskFlow.C");
  AliAnalysisTaskFlowEvent* taskFE = AddTaskFlow(type,METHODS,QA,WEIGHTS);

  
  // Task to check the offline trigger
  if (mode == mLocal || mode == mGRID) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C"); }
  else if (mode == mPROOF || mode == mLocalPAR) {
    gROOT->LoadMacro("AddTaskPhysicsSelection.C"); }
  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();
  if (!DATA) {physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();}
  
 
  // Enable debug printouts
  mgr->SetDebugLevel(2);


  //____________________________________________//
  // Run the analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  if (mode == mLocal || mode == mLocalPAR) {
    mgr->StartAnalysis("local",chain);
  }
  else if (mode == mPROOF) {
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  }
  else if (mode == mGRID) { 
    mgr->StartAnalysis("grid");
  }
  
  timer.Stop();
  timer.Print();
  
}

void CrossCheckUserSettings(Bool_t bData) 
{
 // Check in this method if the user settings make sense.
 
 if(MCEP==kTRUE && bData==kTRUE)
 {
  cout<<endl;
  cout<<"WARNING: In real datasets there is no Monte Carlo information available !!!!"<<endl;
  cout<<"         Set for real data analysis DATA = kTRUE and MCEP = kFALSE in the macro."<<endl;
  cout<<endl;
  exit(0);
 }

} // end of void CrossCheckUserSettings()

void LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal || mode==mGRID) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    cerr<<"libCORRFW loaded..."<<endl;
    if (mode==mLocal) {
      gSystem->Load("libPWG2flowCommon");
      cerr<<"libPWG2flowCommon loaded..."<<endl;
      gSystem->Load("libPWG2flowTasks");
      cerr<<"libPWG2flowTasks loaded..."<<endl;
    }
    if (mode==mGRID) {
      SetupPar("PWG2flowCommon");
      cerr<<"PWG2flowCommon.par loaded..."<<endl;
      SetupPar("PWG2flowTasks");
      cerr<<"PWG2flowTasks.par loaded..."<<endl;
    }
  }
  
  else if (mode == mLocalPAR) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
        SetupPar("CORRFW");
    SetupPar("PWG2flowCommon");
    cerr<<"PWG2flowCommon.par loaded..."<<endl;
    SetupPar("PWG2flowTasks");
    cerr<<"PWG2flowTasks.par loaded..."<<endl;
  }
  
  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    //
    //gEnv->SetValue("XSec.GSI.DelegProxy","2");    
    //  set to debug root versus if needed
    //TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a_dbg");
    //TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a");
    //TProof::Reset("proof://snelling@alicecaf.cern.ch");     
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    // Put appropriate username here
    //TProof::Open("abilandz@alicecaf.cern.ch");
    //TProof::Open("nkolk@alicecaf.cern.ch");
    //TProof::Open("snelling@localhost");
    TProof::Open("alice-caf.cern.ch");
    //TProof::Open("skaf.saske.sk");
    //TProof::Open("prf000-iep-grid.saske.sk");
    //Info("runSKAF.C","Loading libs on proof (may take while, around 1 min) ...");
    // list the data available
    //gProof->ShowDataSets("/*/*"); 
    //gProof->ShowDataSets("/alice/sim/"); //for MC Data
    //gProof->ShowDataSets("/alice/data/"); //for REAL Data 
 
    // Clear the Packages
    /*
    gProof->ClearPackage("STEERBase.par");
    gProof->ClearPackage("ESD.par");
    gProof->ClearPackage("AOD.par");
    gProof->ClearPackage("ANALYSIS.par");
    gProof->ClearPackage("ANALYSISalice.par");
    gProof->ClearPackage("CORRFW.par");
    
    gProof->ClearPackage("PWG2flowCommon");
    gProof->ClearPackage("PWG2flowTasks");
    */
    // Upload the Packages
    gProof->UploadPackage("STEERBase.par");
    gProof->UploadPackage("ESD.par");    
    gProof->UploadPackage("AOD.par");
       
    gProof->UploadPackage("ANALYSIS.par"); 
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->UploadPackage("CORRFW.par");
    gProof->UploadPackage("PWG2flowCommon.par");
    gProof->UploadPackage("PWG2flowTasks.par");

    // Enable the Packages 
    // The global package
    //gProof->EnablePackage("aliroot_v4-19-05-AN",kTRUE);
    // Or separate
    
    gProof->EnablePackage("STEERBase");
    gProof->EnablePackage("ESD");
    gProof->EnablePackage("AOD");
    
    // Always needed
    gProof->EnablePackage("ANALYSIS");
    gProof->EnablePackage("ANALYSISalice");
    gProof->EnablePackage("CORRFW");
    gProof->EnablePackage("PWG2flowCommon");
    gProof->EnablePackage("PWG2flowTasks");

    // Show enables Packages
    gProof->ShowEnabledPackages();
  }  
  
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


// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
    {
      printf("%s not found.\n", aDataDir);
      return 0;
    }
  
  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
    {
      TString execDir(gSystem->pwd());
      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
      TList* dirList            = baseDir->GetListOfFiles();
      Int_t nDirs               = dirList->GetEntries();
      gSystem->cd(execDir);
      
      Int_t count = 0;
      
      for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliESDs.root/esdTree");
	  //  cerr<<presentDirName<<endl;
	}
      
    }
  else
    {
      // Open the input stream
      ifstream in;
      in.open(aDataDir);
      
      Int_t count = 0;
      
      // Read the input list of files and add them to the chain
      TString esdfile;
      while(in.good()) {
	in >> esdfile;
	if (!esdfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
	// add esd file
	chain->Add(esdfile);
      }
      
      in.close();
    }
  
  return chain;
}


// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliAOD.root
  // <aDataDir>/<dir1>/AliAOD.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
    {
      printf("%s not found.\n", aDataDir);
      return 0;
    }
  
  TChain* chain = new TChain("aodTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
    {
      TString execDir(gSystem->pwd());
      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
      TList* dirList            = baseDir->GetListOfFiles();
      Int_t nDirs               = dirList->GetEntries();
      gSystem->cd(execDir);
      
      Int_t count = 0;
      
      for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliAOD.root/aodTree");
	  // cerr<<presentDirName<<endl;
	}
      
    }
  else
    {
      // Open the input stream
      ifstream in;
      in.open(aDataDir);
      
      Int_t count = 0;
      
      // Read the input list of files and add them to the chain
      TString aodfile;
      while(in.good()) {
	in >> aodfile;
	if (!aodfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
	// add aod file
	chain->Add(aodfile);
      }
      
      in.close();
    }
  
  return chain;
}

