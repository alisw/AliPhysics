enum anaModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries
//       (Remark: When using this mode set also Bool_t bUseParFiles = kFALSE; in CreateAlienHandler.C)
//mGridPAR: Analyze files on Grid via AliEn plug-in and using par files for FLOW package
//          (Remark: when using this mode set also Bool_t bUseParFiles = kTRUE; in CreateAlienHandler.C)
 
// CENTRALITY DEFINITION
//Int_t binfirst = 4;  //where do we start numbering bins
//Int_t binlast = 6;  //where do we stop numbering bins
//const Int_t numberOfCentralityBins = 9;
Int_t binfirst = 0;  //where do we start numbering bins
Int_t binlast = 8;  //where do we stop numbering bins
const Int_t numberOfCentralityBins = 9;
Float_t centralityArray[numberOfCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; // in centrality percentile
//Int_t centralityArray[numberOfCentralityBins+1] = {41,80,146,245,384,576,835,1203,1471,10000}; // in terms of TPC only reference multiplicity

TString commonOutputFileName = "outputCentrality"; // e.g.: result for centrality bin 0 will be in the file "outputCentrality0.root", etc


//void runFlowTaskCentralityTrain(Int_t mode=mLocal, Int_t nRuns = 10, 
//Bool_t DATA = kFALSE, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0)

void runFlowTaskCentralityTrain4Candidates(Int_t mode = mGridPAR, Int_t nRuns = 50000000, 
		 Bool_t DATA = kTRUE, const Char_t* dataDir="/alice/data/LHC10h_000137161_p1_plusplusplus", Int_t offset=0) 
//void runFlowTaskCentralityTrain4Candidates(Int_t mode = mLocal, Int_t nRuns = 50000000, 
//				Bool_t DATA = kTRUE, const Char_t* dataDir="./data/", Int_t offset=0) 
//void runFlowTaskCentralityTrain(Int_t mode = mGridPAR, Bool_t DATA = kTRUE)
{
  // Time:
  TStopwatch timer;
  timer.Start();
  // Cross-check user settings before starting:
  //  CrossCheckUserSettings(DATA);
  // Load needed libraries:
  LoadLibraries(mode);
  // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR) 
    {    
      gROOT->LoadMacro("CreateAlienHandler.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandler();  
      if(!alienHandler) return;
    }
  // Chains: 
  if(mode == mLocal || mode == mLocalPAR) {
//    TChain* chain = CreateESDChain(dataDir, nRuns, offset);
    TChain* chain = new TChain();
    chain->Add("~/alice/datasets/Pb/AliESDs.root/esdTree");
    //TChain* chain = CreateAODChain(dataDir, nRuns, offset);
  }
  
  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager"); 
  // Connect plug-in to the analysis manager:
  if(mode == mGrid || mode == mGridPAR) 
    { 
      mgr->SetGridHandler(alienHandler);
    }
  
  // Event handlers:
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  if (!DATA) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc); 
  }

  // Task to check the offline trigger:
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
  AddTaskPhysicsSelection(!DATA);

  //Add the centrality determination task
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AddTaskCentrality();

  // Setup analysis per centrality bin:
  gROOT->LoadMacro("AddTaskFlowCentrality4Candidates.C");
  for (Int_t i=binfirst; i<binlast+1; i++)
  {
    Float_t lowCentralityBinEdge = centralityArray[i];
    Float_t highCentralityBinEdge = centralityArray[i+1];
    Printf("\nWagon for centrality bin %i: %.0f-%.0f",i,lowCentralityBinEdge,highCentralityBinEdge);
    AddTaskFlowCentrality4Candidates( lowCentralityBinEdge,
                                      highCentralityBinEdge,
                                      commonOutputFileName,
                                      0.474, 0.490 );
    AddTaskFlowCentrality4Candidates( lowCentralityBinEdge,
                                      highCentralityBinEdge,
                                      commonOutputFileName,
                                      0.490, 0.506 );
    AddTaskFlowCentrality4Candidates( lowCentralityBinEdge,
                                      highCentralityBinEdge,
                                      commonOutputFileName,
                                      0.506, 0.522 );
  } // end of for (Int_t i=0; i<numberOfCentralityBins; i++)

  // Enable debug printouts:
  mgr->SetDebugLevel(2);
  // Run the analysis:
  if(!mgr->InitAnalysis()) return;  
  mgr->PrintStatus();
  if(mode == mLocal || mode == mLocalPAR) {
    mgr->StartAnalysis("local",chain);
  } else if(mode == mPROOF) {
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  } else if(mode == mGrid || mode == mGridPAR) { 
    mgr->StartAnalysis("grid");
  }

  // Print real and CPU time used for analysis:
  timer.Stop();
  timer.Print();  
  
} // end of void runFlowTaskCentralityTrain(...)

//===============================================================================================
/*
void CrossCheckUserSettings(Bool_t bData) 
{
 // Check in this method if the user settings make sense. 
 if(LYZ1SUM && LYZ2SUM) {cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1 !!!!"<<endl; exit(0); }
 if(LYZ1PROD && LYZ2PROD) {cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1 !!!!"<<endl; exit(0); }
 if(LYZ2SUM && LYZEP) {cout<<" WARNING: You cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2 !!!!"<<endl; exit(0); }
 if(LYZ1SUM && LYZEP) {cout<<" WARNING: You cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2 !!!!"<<endl; exit(0); }
} // end of void CrossCheckUserSettings()
*/
//===============================================================================================

void LoadLibraries(const anaModes mode) 
{
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------

  gSystem->Load("libCore");  
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  
  if (mode==mLocal || mode==mGrid || mode == mGridPAR || mode == mLocalPAR )
  {
    gSystem->Load("libSTEERBase");
    gSystem->Load("libCDB");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libRAWDatarec");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libSTEER");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");   
    gSystem->Load("libTOFbase"); 
    gSystem->Load("libTOFrec"); 

    if (mode == mLocal || mode == mGrid)
    {
      gSystem->Load("libPWGflowBase"); 
      gSystem->Load("libPWGflowTasks"); 
    }
    if (mode == mLocalPAR || mode == mGridPAR )
    {
      AliAnalysisAlien::SetupPar("PWGflowBase");
      AliAnalysisAlien::SetupPar("PWGflowTasks");
    }
  }
  
  //---------------------------------------------------------
  // <<<<<<<<<< PROOF mode >>>>>>>>>>>>
  //---------------------------------------------------------
  else if (mode==mPROOF) {
    //  set to debug root versus if needed
    //TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a_dbg");
    //TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a");
    //TProof::Reset("proof://snelling@alicecaf.cern.ch");     
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof::Open("mkrzewic@alice-caf.cern.ch");
    //TProof::Open("mkrzewic@skaf.saske.sk");
     // list the data available
    //gProof->ShowDataSets("/*/*"); 
    //gProof->ShowDataSets("/alice/sim/"); //for MC Data
    //gProof->ShowDataSets("/alice/data/"); //for REAL Data 
 
    // Clear the Packages
    /*    
    gProof->ClearPackage("STEERBase.par");
    gProof->ClearPackage("ESD.par");
    gProof->ClearPackage("AOD.par");
    */
    //gProof->ClearPackage("ANALYSIS.par");
    //gProof->ClearPackage("ANALYSISalice.par");
    //gProof->ClearPackage("CORRFW.par");
    
    gProof->ClearPackage("PWGflowBase");
    gProof->ClearPackage("PWGflowTasks");
    
    // Upload the Packages
    //gProof->UploadPackage("STEERBase.par");
    //gProof->UploadPackage("ESD.par");    
    //gProof->UploadPackage("AOD.par");
       
    //gProof->UploadPackage("ANALYSIS.par"); 
    //gProof->UploadPackage("ANALYSISalice.par");
    gProof->UploadPackage("CORRFW.par");
    gProof->UploadPackage("PWGflowBase.par");
    gProof->UploadPackage("PWGflowTasks.par");
    gProof->UploadPackage("ALIRECO.par");

    // Enable the Packages 
    // The global package
    TList* list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES","RAW:OCDB:STEER:TOF"));
    gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-07-AN",list);
    gProof->EnablePackage("ALIRECO");
    //gProof->EnablePackage("ANALYSIS");
    //gProof->EnablePackage("ANALYSISalice");
    //gProof->EnablePackage("CORRFW");
    gProof->EnablePackage("PWGflowBase");
    gProof->EnablePackage("PWGflowTasks");

    // Show enables Packages
    gProof->ShowEnabledPackages();
  }  
  
} // end of void LoadLibraries(const anaModes mode) 

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

