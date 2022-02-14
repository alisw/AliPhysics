enum anaModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries
//mGridPAR: Analyze files on Grid via AliEn plug-in and using par files for FLOW package

//const TString analysisType = "AOD"; //"MC", "ESD", "AOD"
const TString analysisType = "MC"; //"MC", "ESD", "AOD"

void runTaskLcToD0MC(Int_t iGroup = 1,
		     TString lDirectory = "/dcache/alice/panosch/alice/sim/2011/LHC11a/CR/Set",
		     Int_t mode = mLocal, Bool_t DATA = kTRUE) {
  //void runSphericityTask(Int_t mode = mGrid, Bool_t DATA = kTRUE) {
  // Time:
  TStopwatch timer;
  timer.Start();
  
  // Load needed libraries:
  LoadLibraries(mode);
  
  // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler();  
    if (!alienHandler) return;
    gROOT->LoadMacro("AliAnalysisTaskSphericityAOD.cxx++");
  }
  // Chains:   
  if(mode==mLocal || mode == mLocalPAR) {
    if (analysisType=="ESD") { 
      TChain* chain = new TChain("esdTree");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set1/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set2/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set3/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set4/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set5/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set6/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set7/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set8/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set9/AliESDs.root");
      chain->Add("/data/alice2/pchrist/HeavyIons/Data/Pass1_4Plus/Set10/AliESDs.root");
    }
    else if(analysisType == "MC") {
      gROOT->LoadMacro("AliAnalysisTaskLcToD0MC.cxx++");
      TChain *chain = new TChain("TE");
      for(Int_t j = 0; j < 10; j++) {
	chain->Add(Form("%s%d/galice.root",lDirectory.Data(),iGroup*(j+1)));
      }
      
      /*for(Int_t i = 0; i < 10; i++) {
	for(Int_t k = 0; k < 10; k++) {
	  for(Int_t l = 0; l < 10; l++) {
	    for(Int_t m = 0; m < 10; m++) { 
	      chain->Add(Form("/dcache/alice/panosch/alice/sim/2011/LHC11h/244340/%d%d%d%d/galice.root",i,k,l,m));
	    }
	  }
	}
	}*/
    }
    else if (analysisType=="AOD") { 
      gROOT->LoadMacro("AliAnalysisTaskSphericityAOD.cxx++");
      TChain* chain = new TChain("aodTree");
      for(Int_t i = 0; i < 10; i++) 
	chain->Add(Form("/dcache/alice/panosch/alice/data/2015/LHC15n/000244340/010%d/AliAOD.root",i));
    }
    else {
      Printf("Selected analysis type not supported!!! Aborting...");
      abort();
    }
  }//local analysis
  //Proof
  if(mode == mPROOF) {
    gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
    gProof->Load("AliAnalysisTaskSphericityAOD.cxx++");
  }

  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("TaskManager");
  // Connect plug-in to the analysis manager:
  if(mode == mGrid || mode == mGridPAR) { 
    mgr->SetGridHandler(alienHandler);
  }
  // Event handlers:
  if(analysisType == "ESD") {
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
  } // end of if(analysisType == "ESD")  
  if(analysisType == "AOD") {
    AliVEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH); 
  } // end of if(analysisType == "AOD")  
  if(analysisType == "MC") {
    AliVEventHandler* esdH = new AliESDInputHandler;
    mgr->SetInputEventHandler(esdH);
    AliMCEventHandler *mc = new AliMCEventHandler();
    mc->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mc); 
  } // end of  if(analysisType == "MC")
    
  // Task to check the offline trigger:
  //if(mode == mLocal || mode == mGrid || mode == mGridPAR)
  if(analysisType == "ESD") {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C"); 
    AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(!DATA);
    if(!DATA){physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();}
    // Enable debug printouts:
    mgr->SetDebugLevel(2);
    
    //Add the centrality determination task
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    //taskCentrality->SelectCollisionCandidates(AliVEvent::kMB);
  }

  // Load the analysis task:
  if(analysisType == "AOD") {
    gROOT->LoadMacro("AddTaskSphericityAOD.C");
    AliAnalysisTaskSphericityAOD* task = AddTaskSphericityAOD();
  }
  else if(analysisType == "MC") {
    gROOT->LoadMacro("AddTaskLcToD0MC.C");
    AliAnalysisTaskLcToD0MC* task = AddTaskLcToD0MC();
  }

  // Run the analysis:
  if(!mgr->InitAnalysis()){return;}
  mgr->PrintStatus(); 
  if(mode == mLocal || mode == mLocalPAR) 
    mgr->StartAnalysis("local",chain);
  else if(mode == mPROOF) 
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  else if(mode == mGrid || mode == mGridPAR) 
    mgr->StartAnalysis("grid");

  // Print real and CPU time used for analysis:  
  timer.Stop();
  timer.Print();
  
} // end of void runSphericityTask(...)

//=============================================================//
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
  if (mode==mLocal || mode==mGrid || mode == mGridPAR) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libEventMixing.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWGTools.so");
    // Use AliRoot includes to compile our task                                   
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    if(mode==mLocal || mode==mGrid)
      gSystem->Load("libPWGCFebye");
    if(mode==mGridPAR)
      SetupPar("PWGCFebye");
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
    //SetupPar("PWGCFebye");
}
  
  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    // Put appropriate username here
    TProof::Open("alice-caf.cern.ch");
    //TProof::Open("skaf.saske.sk");
    //TProof::Open("prf000-iep-grid.saske.sk");

    gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-12-AN");
 
    // Clear the Packages    
    //gProof->ClearPackage("STEERBase.par");
    //gProof->ClearPackage("ESD.par");
    //gProof->ClearPackage("AOD.par");
    //gProof->ClearPackage("ANALYSIS.par");
    //gProof->ClearPackage("ANALYSISalice.par");    
    //gProof->ClearPackage("PWG2ebye");
    
    // Upload the Packages
    //gProof->UploadPackage("STEERBase.par");
    //gProof->UploadPackage("ESD.par");    
    //gProof->UploadPackage("AOD.par");       
    //gProof->UploadPackage("ANALYSIS.par"); 
    //gProof->UploadPackage("ANALYSISalice.par");
    //gProof->UploadPackage("CORRFW.par");
    //gProof->UploadPackage("PWG2ebye");

    // Enable the Packages 
    //gProof->EnablePackage("STEERBase");
    //gProof->EnablePackage("ESD");
    //gProof->EnablePackage("AOD");
    //gProof->EnablePackage("ANALYSIS");
    //gProof->EnablePackage("ANALYSISalice");
    //gProof->EnablePackage("PWG2ebye");

    // Show enables Packages
    //gProof->ShowEnabledPackages();
  }  
  
} // end of void LoadLibraries(const anaModes mode)

//===============================================================================================

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

} // end of void SetupPar(char* pararchivename) 

//===============================================================================================

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

} // end of TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)

//===============================================================================================

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

} // end of TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)

