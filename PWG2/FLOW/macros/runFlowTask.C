enum anaModes {mLocal,mLocalPAR,mPROOF,mGRID};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF

// RUN SETTINGS

// Flow analysis method can be:(set to kTRUE or kFALSE)
Bool_t SP     = kTRUE;
Bool_t LYZ1   = kTRUE;
Bool_t LYZ2   = kFALSE;
Bool_t LYZEP  = kFALSE;
Bool_t GFC    = kTRUE;
Bool_t QC     = kTRUE;
Bool_t FQD    = kTRUE;
Bool_t MCEP   = kTRUE; //not for pp 

Bool_t METHODS[] = {SP,LYZ1,LYZ2,LYZEP,GFC,QC,FQD,MCEP};

// Analysis type can be ESD, AOD, MC, ESDMC0, ESDMC1
const TString type = "ESD";

// Boolean to fill/not fill the QA histograms
Bool_t QA = kTRUE;   

// Boolean to use/not use weights for the Q vector
Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)


//void runFlowTask(Int_t mode=mLocal, Int_t nRuns = 1, 
//	     const Char_t* dataDir="/data/alice2/kolk/PP/LHC09a4/81119", Int_t offset = 0)
void runFlowTask(Int_t mode=mPROOF, Int_t nRuns = 10000, 
		 //	 const Char_t* dataDir="/COMMON/COMMON/LHC09a4_run8158X", Int_t offset = 0)
		 const Char_t* dataDir="/PWG2/akisiel/Therminator_c2030", Int_t offset=0)
{
  TStopwatch timer;
  timer.Start();
  
  LoadLibraries(mode);
  
  if (mode==mLocal || mode == mLocalPAR || mode == mGRID) {
    if (type!="AOD") { TChain* chain = CreateESDChain(dataDir, nRuns, offset);}
    else { TChain* chain = CreateAODChain(dataDir, nRuns, offset);}
  }
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  
  if (type == "ESD"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
  }
  
  if (type == "AOD"){
    AliVEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH); 
  }
  
  if (type == "MC" || type == "ESDMC0" || type == "ESDMC1"){
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
    
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc); }
  
  
  //____________________________________________//
  // Load the tasks
  gROOT->LoadMacro("AddTaskFlow.C");
  AliAnalysisTaskFlowEvent* taskFE = AddTaskFlow(type,METHODS,QA,WEIGHTS);
    
  //____________________________________________//
  // Run the analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  if (mode==mLocal || mode == mLocalPAR) {
    mgr->StartAnalysis("local",chain);
  }
  else if (mode==mPROOF) {
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  }
  else if (mode==mGRID) { 
    mgr->StartAnalysis("local",chain);
  }
  
  timer.Stop();
  timer.Print();
  
}


void LoadLibraries(const anaModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW.so");
    cerr<<"libCORRFW.so loaded..."<<endl;
    gSystem->Load("libPWG2flowCommon.so");
    cerr<<"libPWG2flowCommon.so loaded..."<<endl;
    gSystem->Load("libPWG2flowTasks.so");
    cerr<<"libPWG2flowTasks.so loaded..."<<endl;
  }
  
  else if (mode == mLocalPAR || mode == mGRID) {
    //--------------------------------------------------------
    //If you want to use root and par files from aliroot
    //--------------------------------------------------------  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG2AOD");
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
    
    //  set to debug root versus if needed
    //  TProof::Mgr("alicecaf")->SetROOTVersion("v5-21-01-alice_dbg");
    //  TProof::Mgr("alicecaf")->SetROOTVersion("v5-21-01-alice");
    
    // Connect to proof
    // Put appropriate username here
    // TProof::Reset("proof://snelling@alicecaf.cern.ch"); 
    printf("*** Connect to PROOF ***\n");
    //  TProof::Open("abilandz@alicecaf.cern.ch");
    //    TProof::Open("nkolk@alicecaf.cern.ch");
    TProof::Open("snelling@localhost");
    
    // Enable the STEERBase Package
    //    gProof->ClearPackage("STEERBase.par");
    gProof->UploadPackage("STEERBase.par");
    gProof->EnablePackage("STEERBase");
    // Enable the ESD Package
    //    gProof->ClearPackage("ESD.par");
    gProof->UploadPackage("ESD.par");
    gProof->EnablePackage("ESD");
    // Enable the AOD Package
    //    gProof->ClearPackage("AOD.par");
    gProof->UploadPackage("AOD.par");
    gProof->EnablePackage("AOD");
    // Enable the Analysis Package
    //    gProof->ClearPackage("ANALYSIS.par");
    gProof->UploadPackage("ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
    // Enable the Analysis Package alice
    //    gProof->ClearPackage("ANALYSISalice.par");
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    // Load the PWG2 AOD
    //    gProof->ClearPackage("PWG2AOD.par");
    gProof->UploadPackage("PWG2AOD.par");
    gProof->EnablePackage("PWG2AOD");
    // Enable the Correction Framework
    //    gProof->ClearPackage("CORRFW.par");
    gProof->UploadPackage("CORRFW.par");
    gProof->EnablePackage("CORRFW");
    // Enable Flow Analysis
    //    gProof->ClearPackage("PWG2flowCommon");
    gProof->UploadPackage("PWG2flowCommon.par");
    gProof->EnablePackage("PWG2flowCommon");
    //    gProof->ClearPackage("PWG2flowTasks");
    gProof->UploadPackage("PWG2flowTasks.par");
    gProof->EnablePackage("PWG2flowTasks");
    //
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

