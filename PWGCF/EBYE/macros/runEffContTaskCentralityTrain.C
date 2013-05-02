enum anaModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries
//mGridPAR: Analyze files on Grid via AliEn plug-in and using par files for FLOW package

//Analysis modes
const TString analysisMode = "TPC"; //"TPC", "Global"

//Centrality stuff
Int_t binfirst = 0;  //where do we start numbering bins
Int_t binlast = 8;  //where do we stop numbering bins
const Int_t numberOfCentralityBins = 9;
//Float_t centralityArray[numberOfCentralityBins+1] = {0.,100.}; // in centrality percentile
Float_t centralityArray[numberOfCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; // in centrality percentile
const TString  centralityEstimator = "V0M";
Double_t vertexZ = 10.;
Int_t AODfilterBit = 128;

//output file
TString commonOutputFileName = "AnalysisResults";

//void runEffContTaskCentralityTrain(Int_t mode = mPROOF,Int_t nRuns = 600000, Bool_t DATA = kFALSE, 
                                    //const Char_t* dataDir="/alice/sim/LHC11a10a_000138662_AOD048", Int_t offset=0) { //PbPbAOD
				    //const Char_t* dataDir="/alice/sim/LHC11a10a_000139510", Int_t offset=0) {	    //PbPb ESD			      
void runEffContTaskCentralityTrain(Int_t mode = mLocal, Bool_t DATA = kFALSE) {
  //void runEffContTaskCentralityTrain(const char* runListFileName = "listOfRuns.txt",Int_t mode = mGrid,  Bool_t DATA = kFALSE) {`
  // Time:
  TStopwatch timer;
  timer.Start();
  
  // Load needed libraries:
  LoadLibraries(mode);
  
  // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(runListFileName);  
    if (!alienHandler) return;
    gROOT->LoadMacro("AliAnalysisTaskEffContBF.cxx++");
  }
  // Chains:   
  if(mode==mLocal || mode == mLocalPAR) {
    gROOT->LoadMacro("AliAnalysisTaskEffContBF.cxx++");
    //TChain* chain = new TChain("esdTree");
    //chain->Add("/project/alice/users/alisrm/AcceptanceFilter/MonteCarlo/Set1/AliESDs.root");
    TChain* chain = new TChain("aodTree");
    for (Int_t i=0;i<99; i++){
      filename = "/project/alice/users/alisrm/Efficiency_Contamination/LHC13b3_HIJING_pA_AOD/";
      filename += i; 
      filename += "/AliAOD.root";
      chain->Add(filename.Data());
    }
    //chain->Add("/project/alice/users/alisrm/Efficiency_Contamination/AOD_LHC11a10a/AliAOD.root"); //si da
    //chain->Add("/project/alice/users/alisrm/Efficiency_Contamination/AOD_LHC10d2/AliAOD.root"); 
    //chain->Add("/project/alice/users/pchrist/Data/2011/Set3/AliAOD.root"); //data
  }
  //Proof
  if(mode == mPROOF) {
    gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
    gProof->Load("AliAnalysisTaskEffContBF.cxx++");
  }

  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("FluctuationsAnalysisManager");
  // Connect plug-in to the analysis manager:
  if(mode == mGrid || mode == mGridPAR) { 
    mgr->SetGridHandler(alienHandler);
  }
  
  AliVEventHandler* aodH = new AliAODInputHandler; //NUEVO
  mgr->SetInputEventHandler(aodH); //NUEVO

  //AliMCEventHandler *mc = new AliMCEventHandler();
  //mc->SetReadTR(kFALSE);
  //mgr->SetMCtruthEventHandler(mc); 
    
  // Task to check the offline trigger:
  //if(mode == mLocal || mode == mGrid || mode == mGridPAR)
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C"); 
  //AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(!DATA);
  //physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();
  // Enable debug printouts:
  //mgr->SetDebugLevel(2);
  
  //Add the centrality determination task
  /*gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
  centralityTask->SetMCInput();*/ //antes for ESD

  //centralityTask->SetPass(2);
  //AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  //taskCentrality->SelectCollisionCandidates(AliVEvent::kMB);

  // Load the analysis task:
  gROOT->LoadMacro("AddTaskBalanceEffCont.C");
  //gROOT->LoadMacro("AddEfficiencyTaskCentralityTrain.C");

  //for (Int_t i=binfirst; i<binlast+1; i++) {
  for (Int_t i = 0; i < numberOfCentralityBins; i++) {
    Float_t lowCentralityBinEdge = centralityArray[i];
    Float_t highCentralityBinEdge = centralityArray[i+1];
    Printf("\nWagon for centrality bin %i: %.0f-%.0f",i,lowCentralityBinEdge,highCentralityBinEdge);
    // AddEfficiencyTaskCentralityTrain(analysisMode.Data(),lowCentralityBinEdge, highCentralityBinEdge,commonOutputFileName);
    // AddTaskBalanceEfficiency(kUseHybrid,centralityEstimator, lowCentralityBinEdge, highCentralityBinEdge, vertexZ, AODfilterBit, commonOutputFileName);
  AddTaskBalanceEffCont(centralityEstimator, lowCentralityBinEdge, highCentralityBinEdge, vertexZ, AODfilterBit, commonOutputFileName);  
  } // end of for (Int_t i=0; i<numberOfCentralityBins; i++)

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
  
} // end of void runTaskFluctuations(...)

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
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGTools");
    
    // Use AliRoot includes to compile our task                                   
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
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

    // gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-06-AN");
    gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-44-AN");

    TString extraLibs = "";
    extraLibs += "CORRFW:PWGTools";

    TList *list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS",extraLibs.Data()));

    //gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-06-AN",list);
    gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-44-AN");
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

