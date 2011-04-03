enum anaModes {mLocal,mPROOF,mGrid};
//mLocal: Analyze locally files in your computer using aliroot
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries

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

//void runFlowTaskCentralityPIDTrain(Int_t mode=mLocal, Int_t nEvents = 10,
//Bool_t DATA = kFALSE, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0)

void runFlowTaskCentralityPIDTrain( Int_t mode = mGrid,
                                    Bool_t useFlowParFiles = kTRUE,
                                    Bool_t DATA = kTRUE,
                                    //const Char_t* dataDir="/data/alice2/ab/grid/21-16/TEST/data/",
                                    const Char_t* dataDir="fileList",
                                    Int_t nEvents = 1e4,
                                    Int_t offset=0 )
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

  // Chains:
  if(mode == mLocal)
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
    TChain* chain = CreateESDChain(dataDir, nEvents, offset);
    //TChain* chain = CreateAODChain(dataDir, nEvents, offset);
  }

  // Connect plug-in to the analysis manager:
  if(mode == mGrid)
  {
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(useFlowParFiles);
    if(!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  // Event handlers:
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  if (!DATA)
  {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }

  // Task to check the offline trigger:
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(!DATA);

  //Add the centrality determination task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask* centSelTask = AddTaskCentrality();
  if (!DATA) centSelTask->SetMCInput();


  //Add the TOF tender
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FLOW/macros/AddTaskTenderTOF.C");
  AddTaskTenderTOF();

  // Setup analysis per centrality bin:
  gROOT->LoadMacro("AddTaskFlowCentralityPID.C");
  for (Int_t i=binfirst; i<binlast+1; i++)
  {
    Float_t lowCentralityBinEdge = centralityArray[i];
    Float_t highCentralityBinEdge = centralityArray[i+1];
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kPion,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kPion,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kKaon,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kKaon,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kProton,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kProton,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,2,kTRUE );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kPion,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,3 );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kPion,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,3 );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kKaon,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,3 );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kKaon,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,3 );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kProton,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              -1,3 );
    AddTaskFlowCentralityPID( lowCentralityBinEdge,
                              highCentralityBinEdge,
                              commonOutputFileName,
                              AliPID::kProton,
                              AliFlowTrackCuts::kTOFbetaSimple,
                              1,3 );
  } // end of for (Int_t i=0; i<numberOfCentralityBins; i++)

  // Enable debug printouts:
  mgr->SetDebugLevel(2);
  // Run the analysis:
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(mode == mLocal)
  {
    mgr->StartAnalysis("local",chain);
  }
  else if(mode == mPROOF)
  {
    mgr->StartAnalysis("proof",dataDir,nEvents,offset);
  }
  else if(mode == mGrid)
  {
    mgr->StartAnalysis("grid");
  }

  // Print real and CPU time used for analysis:
  timer.Stop();
  timer.Print();

} // end of void runFlowTaskCentralityPIDTrain(...)

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
    gSystem->Load("libTOFbase");
    gSystem->Load("libTOFrec");
    gSystem->Load("libTRDbase");
    gSystem->Load("libVZERObase");
    gSystem->Load("libVZEROrec");
    gSystem->Load("libT0base");
    gSystem->Load("libT0rec");
    gSystem->Load("libTENDER");
    gSystem->Load("libTENDERSupplies");

    if (useFlowParFiles)
    {
      AliAnalysisAlien::SetupPar("PWG2flowCommon");
      AliAnalysisAlien::SetupPar("PWG2flowTasks");
    }
    else
    {
      gSystem->Load("libPWG2flowCommon");
      gSystem->Load("libPWG2flowTasks");
    }
  }
  else if (mode==mPROOF)
  {
    TList* list = new TList();
    list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
    if (useFlowParFiles)
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:TENDER:TENDERSupplies"));
    else
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:TENDER:TENDERSupplies:PWG2flowCommon:PWG2flowTasks"));

    //list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES","PWG2/FLOW/AliFlowCommon:PWG2/FLOW/AliFlowTasks"));

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
      gProof->UploadPackage("PWG2flowCommon.par");
      gProof->UploadPackage("PWG2flowTasks.par");
    }

    // Show enables Packages
    gProof->ShowEnabledPackages();
  }
} // end of void LoadLibraries(const anaModes mode)

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
    while(in.good())
    {
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

