// $Id$

// runEMCalJetAnalysis.C
// =====================
// This macro can be used to run a jet analysis within the EMCal Jet Framework.
//
// Examples:
// -> Analyze ESDs from the pA pilot run on the AliEn grid with your task in AnaClass.cxx/.h
//     dataType = "ESD", useGrid = kTRUE, pattern = "*ESDs/pass2/*ESDs.root", addCXXs = "AnaClass.cxx", 
//     addHs = "AnaClass.h", gridDir = "/alice/data/2012/LHC12g", gridMode = "full", runNumbers = "188359 188362"
//     
// -> Analyze AODs (up to 96 files) locally given in files_aod.txt
//     dataType = "AOD", useGrid = kFALSE, numLocalFiles = 96
//
// MERGING ON ALIEN
// ++++++++++++++++
// If you run on the grid, you can monitor the jobs with alimonitor.cern.ch. When enough of them are in DONE state,
// you have to merge the output. This can be done automatically, if you just change the gridMode to "terminate" and
// give the EXACT name of the task whose output should be merged in uniqueName.
// 
//
// Authors: R. Haake, S. Aiola

#include <ctime>
#include "TGrid.h"

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, const Int_t nrunspermaster, 
                                     const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, 
                                     Int_t workerTTL, Bool_t isMC);
                                    
//______________________________________________________________________________
void AnalysisTrainCorrJetsLocal (
         const char*    dataType            = "AOD",                       // set the analysis type, AOD, ESD or sESD
         Bool_t         useGrid             = kTRUE,                      // local or grid
	 TString localfilename = "file_aodlhc10d.txt",
         const char*    gridMode            = "test",                      // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
         const char*    pattern             = "*ESDs/pass2/AOD099/*AliAOD.root", //"*/*/AliAOD.root" "*ESDs/pass1/AOD106/*AliAOD.root",    // file pattern (here one can specify subdirs like passX etc.) (used on grid)
         const char*    gridDir             = "/alice/data/2010/LHC10d",   // /alice/data/2011/LHC11d /alice/sim/2012/LHC12f2b   dir on alien, where the files live (used on grid)
         const char*    runNumbers          = /*126437*/" 126432 126425 126424 126422 126409 126408 126407 126406 126405 126404 126403 126359 126352 126351 126350 126285 126284 126283 126168 126167 126160 126158 126097 126090 126088 126082 126081 126078 126073 126008 126007 126004 125855 125851 125850 125849 125848 125847 125844 125843 125842 125633 125632 125630 125296 125134 125101 125100 125097 125085 125023 124751 122375 122374",             // considered run numbers (used on grid) /*LHC12g 188359 188362, LHC11a 146860 146859*/ /*LHC12f2b 158285 159582 */
	 const Int_t nrunspermaster= 100,
         UInt_t         numLocalFiles       = 3,                          // number of files analyzed locally  
         const char*    runPeriod           = "LHC10d",                    // set the run period (used on grid)
         const char*    uniqueName          = "DJetNewCode",     // sets base string for the name of the task on the grid
         UInt_t         pSel                = AliVEvent::kAny,             // used event selection for every task except for the analysis tasks
         Bool_t         useTender           = kFALSE,                      // trigger, if tender task should be used
         Bool_t         isMC                = kFALSE,                      // trigger, if MC handler should be used

         // Here you have to specify additional code files you want to use but that are not in aliroot
         const char*    addCXXs             = "AliAnalysisTaskRecoJetCorrelations.cxx AliAnalysisTaskSEDmesonsForJetCorrelations.cxx",
         const char*    addHs               = "AliAnalysisTaskRecoJetCorrelations.h AliAnalysisTaskSEDmesonsForJetCorrelations.h",

         // These two settings depend on the dataset and your quotas on the AliEN services
         Int_t          maxFilesPerWorker   = 20,
         Int_t          workerTTL           = 86000,
	 Int_t          nfiletestmode       = 1
         )
{

  // Some pre-settings and constants

  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};
  gSystem->SetFPEMask();
  gSystem->Setenv("ETRAIN_ROOT", ".");
  gSystem->Setenv("ETRAIN_PERIOD", runPeriod);
  // change this objects to strings
  TString usedData(dataType);
  TString additionalCXXs(addCXXs);
  TString additionalHs(addHs);
  cout << dataType << " analysis chosen" << endl;
  if (useGrid)  
  {
    cout << "-- using AliEn grid.\n";
    if (usedData == "sESD") 
    {
      cout << "Skimmed ESD analysis not available on the grid!" << endl;
      return;
    }
  }
  else
    cout << "-- using local analysis.\n";
  

  // Load necessary libraries
  LoadLibs();

  // Create analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(uniqueName);

  // Check type of input and create handler for it
  TString localFiles("-1");
  if(usedData == "AOD")
  {
    localFiles = localfilename;
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler* aodH = AddAODHandler();
  }
  else if((usedData == "ESD") || (usedData == "sESD"))
  {
    if (usedData == "ESD")
      localFiles = "files_esd.txt";
    else
      localFiles = "files_sesd.txt";
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler* esdH = AddESDHandler();
  }
  else
  {
    cout << "Data type not recognized! You have to specify ESD, AOD, or sESD!\n";
  }
  
  if(!useGrid)
    cout << "Using " << localFiles.Data() << " as input file list.\n";

  gROOT->LoadMacro("AddTasksCorrJets.C");
  //List of arguments:

// const Int_t iCandtype=1 /*0 = D0, 1=Dstar...*/,
// const TString sCutFile = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root",
// const Double_t dJetPtCut   = 1.,
// const Double_t dJetAreaCut = 0.,
// const Int_t iAccCut = 1,
// const TString sRunPeriod = "LHC10b",
// const Int_t    uBeamType   = 0,
// const UInt_t uTriggerMask = AliVEvent::kMB, /*for jets and phys sel; the D mesons trigger is defined in the cut object*/
// const Bool_t bIsMC=kFALSE,
// TString sText=""/*completes the name of the candidate task lists*/

     AddTasksCorrJets(1,"/data/Work/jets/testEMCalJetFramework/CutFilesMB/DStartoKpipiCuts_new.root",10.,0.,1,runPeriod,0,pSel,isMC,"");

  // Set the physics selection for all given tasks
  TObjArray *toptasks = mgr->GetTasks();
  for (Int_t i=0; i<toptasks->GetEntries(); ++i) 
  {
    AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
    if (!task)
      continue;
    if (task->InheritsFrom("AliPhysicsSelectionTask"))
      continue;
    ::Info("setPSel", "Set physics selection for %s (%s)", task->GetName(), task->ClassName());
    task->SelectCollisionCandidates(pSel);
  }

  if(gridMode=="full") mgr->SetUseProgressBar(1, 25);
        
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();

  if (useGrid) 
  {  // GRID CALCULATION

    AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, gridMode, runNumbers, nrunspermaster, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL, isMC);
    plugin->SetNtestFiles(nfiletestmode); 

    mgr->SetGridHandler(plugin);

    // start analysis
    cout << "Starting GRID Analysis...";
    if(gridMode=="test") mgr->SetDebugLevel(10);
    else mgr->SetDebugLevel(0);
    mgr->StartAnalysis("grid");
  }
  else
  {  // LOCAL CALCULATION

    TChain* chain = 0;
    if (usedData == "AOD") 
    {
      Printf("Run Create AOD Chain");
      gROOT->LoadMacro("/data/Work/jets/testEMCalJetFramework/AODchainWithFriend/CreateAODChain.C");
      chain = CreateAODChain(localFiles.Data(), numLocalFiles,0,kTRUE,kTRUE);
      //Printf("Chain Friend has %d entries", ((TTree*)chain->GetFriend())->GetEntriesFast());
    }
    else
    {  // ESD or skimmed ESD
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numLocalFiles);
    }
    
    // start analysis
    cout << "Starting LOCAL Analysis...";
    mgr->SetDebugLevel(10);
    mgr->StartAnalysis("local", chain);
  }
}

//______________________________________________________________________________
void LoadLibs()
{
  // Load common libraries (better too many than too few)
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTOFbase");
  //gSystem->Load("libTOFrec");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libRAWDatarec.so");
  gSystem->Load("libTPCbase.so");
  gSystem->Load("libTPCrec.so");
  gSystem->Load("libITSbase.so");
  gSystem->Load("libITSrec.so");
  gSystem->Load("libTRDbase.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libSTAT.so");
  gSystem->Load("libTRDrec.so");
  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGDQdielectron");
  gSystem->Load("libPWGHFhfe");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libPWGHFvertexingHF");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libEMCALraw");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFCorrelationsBase");
  gSystem->Load("libPWGCFCorrelationsDPhi");

  //load CGAL, Fastjet and SISCone
  //gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libSISConePlugin");

  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGJEEMCALJetTasks");


  // include paths
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
  gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers, const Int_t nrunspermaster,
                                     const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, 
                                     Int_t workerTTL, Bool_t isMC)
{
  TDatime currentTime;
  TString tmpName(uniqueName);
  /*
  // Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
  if(strcmp(gridMode, "terminate"))
  {
    tmpName += "_";
    tmpName += currentTime.GetDate();
    tmpName += "_";
    tmpName += currentTime.GetTime();
  }else tmpName +="_20130412_122423";
  */
  TString tmpAdditionalLibs("");
  tmpAdditionalLibs = Form("libTree.so libVMC.so libGeom.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libCDB.so libRAWDatabase.so libSTEER.so libANALYSISalice.so libCORRFW.so libTOFbase.so libRAWDatabase.so libRAWDatarec.so libTPCbase.so libTPCrec.so libITSbase.so libITSrec.so libTRDbase.so libTENDER.so libSTAT.so libTRDrec.so libHMPIDbase.so libPWGPP.so libPWGHFbase.so libPWGDQdielectron.so libPWGHFhfe.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so libEMCALUtils.so libPHOSUtils.so libPWGCaloTrackCorrBase.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libTENDER.so libTENDERSupplies.so libPWGEMCAL.so libPWGGAEMCALTasks.so libPWGTools.so libPWGCFCorrelationsBase.so libPWGCFCorrelationsDPhi.so  libCGAL.so libJETAN.so libfastjet.so libSISConePlugin.so libFASTJETAN.so libPWGJE.so libPWGmuon.so libPWGJEEMCALJetTasks.so %s %s",additionalCode.Data(),additionalHeaders.Data());


  TString macroName("");
  TString execName("");
  TString jdlName("");
  macroName = Form("%s.C", tmpName.Data());
  execName = Form("%s.sh", tmpName.Data());
  jdlName = Form("%s.jdl", tmpName.Data());

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetRunMode(gridMode);

  // Here you can set the (Ali)ROOT version you want to use
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-05");
  plugin->SetAliROOTVersion("v5-04-38-AN");

  plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
  plugin->SetDataPattern(pattern); //dir structure in run directory
  plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
  if (!isMC)
   plugin->SetRunPrefix("000");

  plugin->AddRunList(runNumbers);
  plugin->SetNrunsPerMaster(nrunspermaster);

  plugin->SetGridWorkingDir(Form("%s",tmpName.Data()));
  plugin->SetGridOutputDir("output"); 

  plugin->SetAnalysisSource(additionalCode.Data());
  plugin->SetAdditionalLibs(tmpAdditionalLibs.Data());
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks  -I$ALICE_ROOT/PWGJE  -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/PWGJE/EMCALJetTasks -g");
  plugin->AddExternalPackage("boost::v1_43_0");
  plugin->AddExternalPackage("cgal::v3.6");
  plugin->AddExternalPackage("fastjet::v2.4.2");

  plugin->SetDefaultOutputs(kTRUE);
  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);

  //plugin->SetMergeExcludes("");
  plugin->SetAnalysisMacro(macroName.Data());
  plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
  plugin->SetExecutable(execName.Data());
  plugin->SetTTL(workerTTL);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(jdlName.Data());
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");

  return plugin;
}
