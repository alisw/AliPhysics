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
	 TString localfilename = "/data/Work/jets/testEMCalJetFramework/ptJdistrAug2nd/setMCtasks/map/inheritDev/files_LHC10f7a.txt",
         const char*    gridMode            = "test",                      // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
         const char*    pattern             = "*AOD136/*AliAOD.root", //"*/*/AliAOD.root" "*ESDs/pass1/AOD106/*AliAOD.root",    // file pattern (here one can specify subdirs like passX etc.) (used on grid)
         const char*    gridDir             = "/alice/sim/LHC10f7a",   // /alice/data/2011/LHC11d /alice/sim/2012/LHC12f2b   dir on alien, where the files live (used on grid)
         const char*    runNumbers          = /*114931 */" 115186 115193 115393 115401 115414 116102 116288 116402 116403 116562 116571 116574 116643 116645 117048 117050 117052 117053 117054 117059 117060 117063 117065 117077 117086 117092 117099 117109 117112 117116 117220 117222 119159 119161 119163 119841 119842 119844 119845 119846 119849 119853 119856 119859 119862 120067 120069 120072 120073 120076 120079 120244 120503 120504 120505 120616 120617 120671 120741 120750 120758 120820 120821 120822 120823 120824 120825 120829 122374 122375 124187 124191 124355 124358 124362 124367 124378 124380 124381 124604 124605 124606 124607 124608 124702 124746 124750 124751 125023 125085 125097 125100 125101 125133 125134 125139 125140 125156 125186 125296 125628 125630 125632 125633 125842 125843 125844 125847 125848 125849 125850 125851 125855 126004 126007 126008 126073 126078 126081 126082 126088 126090 126097 126158 126160 126167 126168 126283 126284 126285 126350 126351 126352 126359 126403 126404 126405 126406 126407 126408 126409 126422 126424 126425 126432 126437 127719 127724 127729 127730 127814 127819 127930 127940 128263 128778 128913 129536 129599 129639 129641 129654 129659 129666 129667 129723 129725 129726 129729 129735 129736 129738 129742 129744 129959 129960 129961 129962 129966 129983 130149 130151 130157 130158 130172 130178 130179 130342 130343 130354 130356 130358 130360 130375 130479 130480 130481 130517 130519 130520 130524 130526 130601 130608 130620 130621 130623 130628 130696 130704 130793 130795 130798 130799 130834 130840 130842 130844 130847 130848",             // considered run numbers (used on grid) /*LHC12g 188359 188362, LHC11a 146860 146859*/ /*LHC12f2b 158285 159582 */
	 const Int_t nrunspermaster= 100,
         UInt_t         numLocalFiles       = 10,                          // number of files analyzed locally  
         const char*    runPeriod           = "LHC10f7a",                    // set the run period (used on grid)
         const char*    uniqueName          = "DJetNewCodeMCSandBchJ",     // sets base string for the name of the task on the grid
         UInt_t         pSel                = AliVEvent::kAny,             // used event selection for every task except for the analysis tasks
         Bool_t         useTender           = kFALSE,                      // trigger, if tender task should be used
         Bool_t         isMC                = kTRUE,                      // trigger, if MC handler should be used
	 Bool_t         isReco                = kTRUE,
	 Bool_t         isMap                = kTRUE,
         // Here you have to specify additional code files you want to use but that are not in aliroot
         const char*    addCXXs             = "",
         const char*    addHs               = "",

         // These two settings depend on the dataset and your quotas on the AliEN services
         Int_t          maxFilesPerWorker   = 20,
         Int_t          workerTTL           = 86000,
	 Int_t          nfiletestmode       = 3
         )
{

  // Some pre-settings and constants
  TStopwatch watch;
  watch.Start();

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

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/FlavourJetTasks/macros/AddTasksFlavourJet.C");
  //List of arguments:

  // const Int_t iCandType = 1 /*0 = D0, 1=Dstar...*/,
  //  const TString sCutFile = "cutsHF/D0toKpiCutsppRecVtxNoPileupRejNoEMCAL.root",
  //  const Double_t dJetPtCut   = 1.,
  //  const Double_t dJetAreaCut = 0.,
  //  const char *acctype = "TPC",
  //  const TString sRunPeriod = "LHC10b",
  //  const Int_t    uBeamType = 0,
  //  const UInt_t uTriggerMask = AliVEvent::kMB, /*for jets; the D mesons trigger is defined in the cut object*/
  //  const Bool_t bIsMC = kFALSE,
  //  const Bool_t bIsReco = kFALSE,
  //  const Bool_t bIsMap = kFALSE,
  //  TString sText=""/*completes the name of the candidate task lists*/

  //AddTasksFlavourJet(1,"/data/Work/jets/testEMCalJetFramework/CutFilesMB/DStartoKpipiCuts_new.root",10.,0.,"TPC",runPeriod,0,pSel,isMC,isReco,isMap,"");
  AddTasksFlavourJet(1,"/data/Work/jets/testEMCalJetFramework/CutFilesMB/DStartoKpipiCuts.root",5.,0.,"TPC",runPeriod,0,pSel,isMC,isReco,isMap,"");

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
  watch.Stop();
  watch.Print();
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
  //gSystem->Load("libTOFbase");
  //gSystem->Load("libTOFrec");
  //gSystem->Load("libRAWDatabase");
  //gSystem->Load("libRAWDatarec");
  //gSystem->Load("libTPCbase");
  //gSystem->Load("libTPCrec");
  //gSystem->Load("libITSbase");
  //gSystem->Load("libITSrec");
  //gSystem->Load("libTRDbase");
  //gSystem->Load("libTender");
  //gSystem->Load("libSTAT");
  //gSystem->Load("libTRDrec");
  //gSystem->Load("libHMPIDbase");
  gSystem->Load("libPWGTools");
  //gSystem->Load("libPWGPP");
  gSystem->Load("libPWGHFbase");
  //gSystem->Load("libPWGDQdielectron");
  //gSystem->Load("libPWGHFhfe");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGHFvertexingHF");
  gSystem->Load("libEMCALUtils");
  //gSystem->Load("libPHOSUtils");
  //gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libEMCALraw");
  gSystem->Load("libEMCALbase");
  //gSystem->Load("libEMCALrec");
  //gSystem->Load("libTRDbase");
  //gSystem->Load("libVZERObase");
  //gSystem->Load("libVZEROrec");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libESDfilter");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  //gSystem->Load("libPWGCFCorrelationsBase");
  //gSystem->Load("libPWGCFCorrelationsDPhi");

  //load CGAL, Fastjet3 
  //gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");
  gSystem->Load("libfastjetcontribfragile");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGGAEMCALTasks");

  gSystem->Load("libPWGJEEMCALJetTasks");
  gSystem->Load("libPWGJEFlavourJetTasks");


  // include paths
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PWG/EMCAL -I$ALICE_ROOT/PWGJE/EMCALJetTasks");
  gSystem->AddIncludePath(" -I$ALICE_ROOT/PWGHF/ -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGJE/FlavourJetTasks");
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
  tmpAdditionalLibs = Form("libTree.so libVMC.so libGeom.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libPhysics.so libSTEERBase.so libESD.so libESDfilter.so libAOD.so libOADB.so libANALYSIS.so libCDB.so libRAWDatabase.so libSTEER.so libANALYSISalice.so libCORRFW.so libPWGTools.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so libEMCALUtils.so libEMCALraw.so libEMCALbase.so libTender.so libPWGEMCAL.so libPWGGAEMCALTasks.so libPWGTools.so libCGAL.so libfastjet.so libsiscone.so libsiscone_spherical.so libfastjetplugins.so libfastjettools.so libfastjetcontribfragile.so libJETAN.so libFASTJETAN.so  libPWGJE.so libPWGJEEMCALJetTasks.so libPWGJEFlavourJetTasks.so %s %s",additionalCode.Data(),additionalHeaders.Data());


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
  plugin->SetROOTVersion(":v5-34-08-6");
  plugin->SetAliROOTVersion("vAN-20141005");

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
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks  -I$ALICE_ROOT/PWGJE  -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/PWGJE/EMCALJetTasks -I$ALICE_ROOT/PWGJE/FlavourJetTasks -g");
  plugin->AddExternalPackage("boost::v1_53_0");
  plugin->AddExternalPackage("cgal::v4.4");
  plugin->AddExternalPackage("fastjet::v3.0.6_1.012");

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
