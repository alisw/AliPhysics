enum anaModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries
//       (Remark: When using this mode set also Bool_t bUseParFiles = kFALSE; in CreateAlienHandler.C)
//mGrid + par files: Analyze files on Grid via AliEn plug-in and using par files for FLOW package.
//                   Simply set Int_t mode = mGrid and Bool_t useFlowParFiles = kTRUE as arguments.

// Shi
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TString.h>
#include <TFile.h>


#include "TROOT.h"
#include "TString.h"

#include <TProofMgr.h>
#include <TProof.h>
#endif 

class AliAnalysisGrid;
class AliAnalysisAlien;

// CENTRALITY DEFINITION
Bool_t kUseCentrality = kTRUE;
//Int_t binfirst = -1; //if kUseCentrality then change accordingly
//Int_t binlast = -1;  //if kUseCentrality then change accordingly
Int_t binfirst = 8; //if kUseCentrality then change accordingly
Int_t binlast = 9;  //if kUseCentrality then change accordingly
const Int_t numberOfCentralityBins = 9;
Float_t centralityArray[numberOfCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; // in centrality percentile
//Int_t centralityArray[numberOfCentralityBins+1] = {41,80,146,245,384,576,835,1203,1471,10000}; // in terms of TPC only reference multiplicity

TString commonOutputFileName = "AnalysisResults"; // e.g.: result for centrality bin 0 will be in the file "outputCentrality0.root", etc

void LoadLibraries(const anaModes mode, Bool_t useFlowParFiles );
TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset);
TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset);

//void runFlowTask(Int_t mode=mLocal, Int_t nRuns = 10, 
//Bool_t DATA = kFALSE, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0)
void runFlowTask(Int_t mode=mLocal, Int_t nRuns = 15, Bool_t useFlowParFiles = kFALSE,
		 Bool_t DATA = kTRUE, const Char_t* dataDir="/home/alidock/AODdata", Int_t offset = 0)
//void runFlowTask(Int_t mode = mGridPAR, Int_t nRuns = 50000000, 
//		 Bool_t DATA = kTRUE, const Char_t* dataDir="/alice/data/LHC10h_000137161_p1_plusplusplus", Int_t offset=0) 
//void runFlowTask(Int_t mode = mLocal, Int_t nRuns = 50000000, 
//		 Bool_t DATA = kTRUE, const Char_t* dataDir="./data/", Int_t offset=0) 
//void runFlowTask(Int_t mode = mGridPAR, Bool_t DATA = kTRUE)
//void runFlowTask(Int_t mode = mGrid,
                 //Bool_t useFlowParFiles = kFALSE,
                 //Bool_t DATA = kTRUE,
                 //Bool_t useTender = kFALSE,
                 //const Char_t* dataDir= "",
                 //Int_t offset=0,
                 //Int_t nRuns = 50000000)
{
  // Time:
  TStopwatch timer;
  timer.Start();
  // Cross-check user settings before starting:
  //  CrossCheckUserSettings(DATA);
  // Load needed libraries:
  LoadLibraries(mode,useFlowParFiles);
  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager"); 
  // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR) 
    {    
	  /*gROOT->LoadMacro("CreateAlienHandler.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandler(useFlowParFiles);  
      if(!alienHandler) return;*/
	  #if !defined (__CINT__) || defined (__CLING__)
        gInterpreter->LoadMacro("CreateAlienHandler.C");
        AliAnalysisGrid *alienHandler = reinterpret_cast<AliAnalysisGrid*>(gInterpreter->ProcessLine(Form("CreateAlienHandler(%d)",useFlowParFiles)));  
      #else
        gROOT->LoadMacro("CreateAlienHandler.C");
        AliAnalysisGrid *alienHandler = CreateAlienHandler(useFlowParFiles);  
      #endif
      if(!alienHandler) return;
      else mgr->SetGridHandler(alienHandler);
    }
  // Chains: 
  //TChain *chain = new TChain("esdTree"); 
  TChain *chain = new TChain("aodTree"); 
  if(mode == mLocal || mode == mLocalPAR) {
    //TChain *chain = new TChain("esdTree");
    //chain->Add("/home/pchrist/ALICE/HeavyIons/Data/137161/Set1/AliESDs.root");
    //TChain* chain = CreateESDChain(dataDir, nRuns, offset);
    //chain = CreateESDChain(dataDir, nRuns, offset);
    //TChain* chain = CreateAODChain(dataDir, nRuns, offset);
    chain = CreateAODChain(dataDir, nRuns, offset);
  }
  
  /*
  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager"); 
  // Connect plug-in to the analysis manager:
  if(mode == mGrid || mode == mGridPAR) 
    { 
      mgr->SetGridHandler(alienHandler);
    }
  */
  // Event handlers:
  //AliVEventHandler* esdH = new AliESDInputHandler;
  //mgr->SetInputEventHandler(esdH);
  AliVEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  if (!DATA) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc); 
  }

  // Task to check the offline trigger:
  #if !defined (__CINT__) || defined (__CLING__)
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
    reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form("AddTaskPhysicsSelection(%d)",!DATA)));

    //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    //Bool_t bkgRej = kTRUE;
    //AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gROOT->ProcessLine(Form("AddTaskPhysicsSelection(%d,%d)", isMC,bkgRej)));
    //if (isMC){
    //  physSelTask->GetPhysicsSelection()->SetAnalyzeMC();
    //}
  
  #else
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
    AddTaskPhysicsSelection(!DATA);
  #endif

  //Add the centrality determination task
  if(kUseCentrality) {
	#if !defined (__CINT__) || defined (__CLING__)
	  cout<<"Loading TaskMultSelection"<<endl;
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AliMultSelectionTask *multSelTask = reinterpret_cast<AliMultSelectionTask*>(gROOT->ProcessLine("AddTaskMultSelection(kFALSE,\"A\")")); 
      multSelTask->SetStoreAllQA(kTRUE);
      //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      //reinterpret_cast<AliCentralitySelectionTask *>(gROOT->ProcessLine(Form("AddTaskCentrality()")));
    #else
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AliMultSelectionTask *multSelTask = AddTaskMultSelection(kFALSE,"A");
      multSelTask->SetStoreAllQA(kTRUE);
      //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      //AddTaskCentrality();
    #endif
  }
  
  //Add the TOF tender
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/FLOW/macros/AddTaskTenderFlow.C");
  //AddTaskTenderFlow();

  // Setup analysis and usage of centrality bins
  /*#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AddTaskFlowQCFineCentBinning.C");
  #else
    gROOT->LoadMacro("AddTaskFlowQCFineCentBinning.C");
  #endif
  */
  // Setup analysis and usage of centrality bins
  #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AddTaskCRC.C");
  #else
    gROOT->LoadMacro("AddTaskCRC.C");
  #endif
  
  //#if !defined (__CINT__) || defined (__CLING__)
    //gInterpreter->LoadMacro("AddTaskAccContForWeights.C");
  //#else
    //gROOT->LoadMacro("AddTaskAccContForWeights.C");
  //#endif
  
  //#if !defined (__CINT__) || defined (__CLING__)
    //gInterpreter->LoadMacro("AddTaskFlow.C");
  //#else
    //gROOT->LoadMacro("AddTaskFlow.C");
  //#endif
  
  Float_t kLowCentralityBin = -1.;
  Float_t kHighCentralityBin = -1;
  if(kUseCentrality) {
    kLowCentralityBin = centralityArray[binfirst];
    kHighCentralityBin = centralityArray[binlast];
  }

                      
  //#if !defined (__CINT__) || defined (__CLING__)
    //static_cast<void>(gInterpreter->ProcessLine(Form("AddTaskFlowQCFineCentBinning()")));
  //#else
    //AddTaskFlowQC();
  //#endif
  
  #if !defined (__CINT__) || defined (__CLING__)
    static_cast<void>(gInterpreter->ProcessLine(Form("AddTaskCRC()")));
  #else
    AddTaskCRC();
  #endif
  
  //#if !defined (__CINT__) || defined (__CLING__)
    //static_cast<void>(gInterpreter->ProcessLine(Form("AddTaskAccContForWeights()")));
  //#else
    //AddTaskAccContForWeights();
  //#endif
  
/*  #if !defined (__CINT__) || defined (__CLING__)
    static_cast<void>(gInterpreter->ProcessLine(Form("AddTaskFlow(%f, %f, \"%s\")", kLowCentralityBin, kHighCentralityBin, commonOutputFileName.Data())));
    //TString AddTaskFlowString = Form("AddTaskFlow(%f, %f, %s)", kLowCentralityBin, kHighCentralityBin, commonOutputFileName.Data());
    //static_cast<void>(gInterpreter->ProcessLine(AddTaskFlowString.Data()));
  #else
    AddTaskFlow(kLowCentralityBin,
	        kHighCentralityBin,
	        commonOutputFileName);
  #endif
*/

  // Enable debug printouts:
  mgr->SetDebugLevel(10);
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
  
} // end of void runFlowTask(...)

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

void LoadLibraries(const anaModes mode, Bool_t useFlowParFiles )
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
  gSystem->Load("libMinuit");

  if (mode==mLocal || mode==mGrid)
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
    gSystem->Load("libTPCbase");
    gSystem->Load("libTOFbase");
    gSystem->Load("libTOFrec");
    gSystem->Load("libTRDbase");
    gSystem->Load("libVZERObase");
    gSystem->Load("libVZEROrec");
    gSystem->Load("libT0base");
    gSystem->Load("libT0rec");
    gSystem->Load("libTender");
    gSystem->Load("libTenderSupplies");

    if (useFlowParFiles)
    {
      AliAnalysisAlien::SetupPar("PWGflowBase");
      AliAnalysisAlien::SetupPar("PWGflowTasks");
    }
    else
    {
      gSystem->Load("libPWGflowBase");
      gSystem->Load("libPWGflowTasks");
    }
  }
  else if (mode==mPROOF)
  {
    TList* list = new TList();
    list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
    if (useFlowParFiles)
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:Tender:TenderSupplies"));
    else
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:Tender:TenderSupplies:PWGflowBase:PWGflowTasks"));

    //list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES","PWG/FLOW/Base:PWG/FLOW/Tasks"));

    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    //TProof* proof = TProof::Open("alice-caf.cern.ch");
    TProof* proof = TProof::Open("skaf.saske.sk");

    // list the data available
    //gProof->ShowDataSets("/*/*");
    //gProof->ShowDataSets("/alice/sim/"); //for MC Data
    //gProof->ShowDataSets("/alice/data/"); //for REAL Data

    proof->ClearPackages();
    proof->EnablePackage("VO_ALICE@AliRoot::v4-21-14-AN",list);

    if (useFlowParFiles)
    {
      gProof->UploadPackage("PWGflowBase.par");
      gProof->UploadPackage("PWGflowTasks.par");
    }

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

